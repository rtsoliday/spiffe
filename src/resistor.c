/*************************************************************************\
 * Copyright (c) 2002 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2002 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: resistor.c
 * purpose: process resistor definitions and implement resistors 
 *
 * Michael Borland, 1992 
 */
#include "spiffe.h"

void process_resistor_definition(
                                 RESISTOR *resistor,
                                 FIELDS *EM_problem,
                                 NAMELIST_TEXT *nl_text)
{
#include "resistor.h"
  long nr;

  if ((nr = resistor->n_resistors) == 0)
    {
      resistor->start = resistor->end = resistor->position = resistor->conductivity = NULL;
      resistor->istart = resistor->iend = resistor->iposition = resistor->direction = NULL;
    }
  resistor->start = trealloc(resistor->start, sizeof(*resistor->start) * (nr + 1));
  resistor->end = trealloc(resistor->end, sizeof(*resistor->end) * (nr + 1));
  resistor->istart = trealloc(resistor->istart, sizeof(*resistor->istart) * (nr + 1));
  resistor->iend = trealloc(resistor->iend, sizeof(*resistor->iend) * (nr + 1));
  resistor->position = trealloc(resistor->position, sizeof(*resistor->position) * (nr + 1));
  resistor->iposition = trealloc(resistor->iposition, sizeof(*resistor->iposition) * (nr + 1));
  resistor->direction = trealloc(resistor->direction, sizeof(*resistor->direction) * (nr + 1));
  resistor->conductivity = trealloc(resistor->conductivity, sizeof(*resistor->conductivity) * (nr + 1));

  start = end = position = 0;
  direction = "z";
  conductivity = sqrt(DBL_MAX);

  process_namelist(&define_resistor, nl_text);
  print_namelist(stdout, &define_resistor);

  if (start >= end)
    bomb("start >= end for resistor ", NULL);
  resistor->start[nr] = start;
  resistor->end[nr] = end;

  if (strncmp("r", direction, 1) == 0)
    {
      /* induced current is Jr */
      if (end > EM_problem->rmax || start < 0)
        bomb("start or end out of range for resistor", NULL);
      if (position < EM_problem->zmin || position > EM_problem->zmax)
        bomb("position out of range for resistor", NULL);
      resistor->istart[nr] = start / EM_problem->dr + 0.5;
      resistor->iend[nr] = end / EM_problem->dr + 0.5;
      resistor->iposition[nr] = (position - EM_problem->zmin) / EM_problem->dz;
      resistor->direction[nr] = 1;
    }
  else if (!direction[0] || strncmp("z", direction, 1) == 0)
    {
      /* induced current is Jz */
      if (position > EM_problem->rmax || position < 0)
        bomb("position out of range for resistor", NULL);
      if (start < EM_problem->zmin || end > EM_problem->zmax)
        bomb("start or end out of range for resistor", NULL);
      resistor->istart[nr] = (start - EM_problem->zmin) / EM_problem->dz + 0.5;
      resistor->iend[nr] = (end - EM_problem->zmin) / EM_problem->dz + 0.5;
      resistor->iposition[nr] = position / EM_problem->dr;
      resistor->direction[nr] = 0;
    }
  else
    bomb("unknown direction specified for resistor", NULL);

  if (conductivity < 0)
    bomb("negative conductivity is meaningless", NULL);
  resistor->conductivity[nr] = conductivity;

  resistor->n_resistors += 1;
}

void add_resistor_currents(FIELDS *EM_problem, RESISTOR *resistor)
{
  long i, ir;
  static long first_call = 1;

  for (ir = 0; ir < resistor->n_resistors; ir++)
    {
      if (resistor->conductivity[ir] > (epsilon_o / EM_problem->dti))
        bomb("conductivity too large for resistor--reduce the conductivity or decrease the time step", NULL);
      if (resistor->direction[ir] == 0)
        {
          /* Jz = sigma*Ez */
          for (i = resistor->istart[ir]; i <= resistor->iend[ir]; i++)
            EM_problem->Jz[i][resistor->iposition[ir]] +=
              EM_problem->Ez[i][resistor->iposition[ir]] * resistor->conductivity[ir];
        }
      else
        {
          /* Jr = sigma*Er */
          for (i = resistor->istart[ir]; i <= resistor->iend[ir]; i++)
            EM_problem->Jr[resistor->iposition[ir]][i] +=
              EM_problem->Er[resistor->iposition[ir]][i] * resistor->conductivity[ir];
        }
    }
  if (first_call)
    {
      printf("resistor check:\n");
      for (ir = 0; ir < resistor->n_resistors; ir++)
        {
          printf("#%-2ld: direction %s, conductivity %e, from grid %ld to %ld at grid %ld\n",
                 ir, (resistor->direction[ir] ? "r" : "z"),
                 resistor->conductivity[ir], resistor->istart[ir], resistor->iend[ir],
                 resistor->iposition[ir]);
        }
      first_call = 0;
    }
}
