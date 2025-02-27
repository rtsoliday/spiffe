/*************************************************************************\
 * Copyright (c) 2002 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2002 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: geometry.c
 * purpose: process geometry definition for EM problem 
 *
 * Michael Borland, 1992 
 */
#include "spiffe.h"
#include "table.h"

void process_antenna_definition(
                                ANTENNA **antenna,
                                long n_antennae,
                                FIELDS *EM_problem,
                                NAMELIST_TEXT *nl_text)
{
#include "antenna.h"
  SDDS_DATASET SDDSin;
  long i;

  start = end = position = 0;
  direction = "z";
  current = frequency = time_offset = phase = 0;
  waveform = NULL;

  *antenna = trealloc(*antenna, sizeof(**antenna) * (n_antennae + 1));

  process_namelist(&define_antenna, nl_text);
  print_namelist(stdout, &define_antenna);

  if (start >= end)
    bomb("start >= end for antenna ", NULL);
  (*antenna)[n_antennae].start = start;
  (*antenna)[n_antennae].end = end;
  (*antenna)[n_antennae].position = position;

  if (strncmp("r", direction, 1) == 0)
    {
      /* current is Jr */
      if (start < 0 || end > EM_problem->rmax)
        bomb("antenna start or end out of range", NULL);
      if (position < EM_problem->zmin || position > EM_problem->zmax)
        bomb("antenna position out of range", NULL);
      (*antenna)[n_antennae].istart = start / EM_problem->dr + 0.5;
      (*antenna)[n_antennae].iend = end / EM_problem->dr + 0.5;
      (*antenna)[n_antennae].iposition = (position - EM_problem->zmin) / EM_problem->dz;
      (*antenna)[n_antennae].direction = 1;
    }
  else if (strncmp("phi", direction, 3) == 0)
    {
      /* current is Jphi */
      if (start < EM_problem->zmin || end > EM_problem->zmax)
        bomb("antenna start or end out of range", NULL);
      if (position < 0 || position > EM_problem->rmax)
        bomb("antenna position out of range", NULL);
      (*antenna)[n_antennae].istart = (start - EM_problem->zmin) / EM_problem->dz;
      (*antenna)[n_antennae].iend = (end - EM_problem->zmin) / EM_problem->dz;
      (*antenna)[n_antennae].iposition = position / EM_problem->dr;
      (*antenna)[n_antennae].direction = 2;
    }
  else if (!direction[0] || strncmp("z", direction, 1) == 0)
    {
      /* current is Jz */
      if (start < EM_problem->zmin || end > EM_problem->zmax)
        bomb("antenna start or end out of range", NULL);
      if (position < 0 || position > EM_problem->rmax)
        bomb("antenna position out of range", NULL);
      (*antenna)[n_antennae].istart = (start - EM_problem->zmin) / EM_problem->dz + 0.5;
      (*antenna)[n_antennae].iend = (end - EM_problem->zmin) / EM_problem->dz + 0.5;
      (*antenna)[n_antennae].iposition = position / EM_problem->dr;
      (*antenna)[n_antennae].direction = 0;
    }

  if (!((*antenna)[n_antennae].current = current))
    printf("warning: zero current antenna will be ignored");
  if (((*antenna)[n_antennae].frequency = frequency) <= 0)
    bomb("antenna frequency <= 0", NULL);
  (*antenna)[n_antennae].time_offset = time_offset;
  (*antenna)[n_antennae].phase = phase;

  if (!waveform)
    bomb("no antenna waveform file given", NULL);
  if (!SDDS_InitializeInput(&SDDSin, waveform))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
  if (SDDS_ReadPage(&SDDSin) <= 0 ||
      ((*antenna)[n_antennae].n_pts = SDDS_CountRowsOfInterest(&SDDSin)) <= 2)
    {
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
      bomb("SDDS waveform file empty or invalid", NULL);
    }
  if (!((*antenna)[n_antennae].t_wf = SDDS_GetColumnInDoubles(&SDDSin, "t")) ||
      !((*antenna)[n_antennae].W_wf = SDDS_GetColumnInDoubles(&SDDSin, "W")))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

  (*antenna)[n_antennae].last_index = 0;
  for (i = 0; i < (*antenna)[n_antennae].n_pts; i++)
    (*antenna)[n_antennae].t_wf[i] += time_offset;

  for (i = 1; i < (*antenna)[n_antennae].n_pts; i++)
    if ((*antenna)[n_antennae].t_wf[i - 1] >= (*antenna)[n_antennae].t_wf[i])
      bomb("waveform data is not ordered by increasing time index", NULL);
  if (!SDDS_Terminate(&SDDSin))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
}

void add_antenna_excitation(FIELDS *EM_problem, ANTENNA *antenna, long n_antennae, double time)
{
  double *t, *W = NULL, J;
  long i, j, n_pts = 0;

  for (j = 0; j < n_antennae; j++)
    {
      if (!(t = antenna->t_wf) || !(W = antenna->W_wf) || (n_pts = antenna->n_pts) <= 1)
        bomb("missing antenna arrays", NULL);

      if (time <= t[0])
        J = W[0];
      else if (time > t[n_pts - 1] || antenna->last_index >= (antenna->n_pts - 1))
        J = W[n_pts - 1];
      else
        {
          if ((i = antenna->last_index) != 0 && time < t[i])
            bomb("antenna time indexing error", NULL);
          while (i < n_pts - 2 && time >= t[i + 1])
            i++;
          antenna->last_index = i;
          J = INTERPOLATE(W[i + 1], W[i], t[i + 1], t[i], time);
        }
      J *= antenna->current * sin(PIx2 * antenna->frequency * time + antenna->phase);
      if (antenna->direction == 0)
        {
          /* Jz */
          for (i = antenna->istart; i <= antenna->iend; i++)
            EM_problem->Jz[i][antenna->iposition] += J;
        }
      else if (antenna->direction == 1)
        {
          /* Jr */
          for (i = antenna->istart; i <= antenna->iend; i++)
            EM_problem->Jr[antenna->iposition][i] += J;
        }
      else if (antenna->direction == 2)
        {
          /* Jphi */
          for (i = antenna->istart; i <= antenna->iend; i++)
            EM_problem->Jphi[i][antenna->iposition] += J;
        }
      antenna++;
    }
}
