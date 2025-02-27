/*************************************************************************\
 * Copyright (c) 2002 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2002 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: load_particles.c
 * purpose: setup and execution of particle loading
 *
 * Michael Borland, 1992 
 */
#include "spiffe.h"
#include "match_string.h"

#include "load_particles.h"

#define DEBUG 0

void perform_particle_loading(
                              FIELDS *EM_problem,
                              BEAM *beam,
                              NAMELIST_TEXT *nl_text)
{
  SDDS_DATASET SDDSin;
  long np_new, i, ip;
  double *q, *z = NULL, *r = NULL, *pz = NULL, *pr = NULL, *pphi = NULL, gamma;
  long *particleID = NULL;
  short hasParticleID = 0;
  filename = NULL;
  sample_interval = 1;
  stiffness = 1;

  process_namelist(&load_particles, nl_text);
  print_namelist(stdout, &load_particles);

  if (!filename || !strlen(filename) || !fexists(filename))
    bomb("No file given or found.", NULL);
  if (sample_interval <= 0)
    bomb("sample_interval <= 0", NULL);

  if (!SDDS_InitializeInput(&SDDSin, filename))
    SDDS_Bomb("unable to open input file for particle loading");

  if (SDDS_CheckColumn(&SDDSin, "z", "m", SDDS_DOUBLE, stderr) != SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "r", "m", SDDS_DOUBLE, stderr) != SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "pz", "m$be$nc", SDDS_DOUBLE, stderr) != SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "pr", "m$be$nc", SDDS_DOUBLE, stderr) != SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "pphi", "m$be$nc", SDDS_DOUBLE, stderr) != SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "q", "C", SDDS_DOUBLE, stderr) != SDDS_CHECK_OK)
    {
      SDDS_Bomb("Invalid particle input");
    }
  if (SDDS_CheckColumn(&SDDSin, "particleID", NULL, SDDS_LONG, NULL) == SDDS_CHECK_OK)
    hasParticleID = 1;

  beam->stiffness = stiffness;
  beam->QoMC = -e_mks / (me_mks * c_mks);
  while (SDDS_ReadPageSparse(&SDDSin, 0, sample_interval, 0, 0) > 0)
    {
      if ((np_new = SDDS_RowCount(&SDDSin)) <= 0)
        continue;
#if DEBUG
      fprintf(stderr, "%ld particles in memory now\nread %ld rows from file\n", beam->np, np_new);
#endif
      if (!(beam->Q = SDDS_Realloc(beam->Q, sizeof(*beam->Q) * (beam->np + np_new))) ||
          !(beam->z = SDDS_Realloc(beam->z, sizeof(*beam->z) * (beam->np + np_new))) ||
          !(beam->r = SDDS_Realloc(beam->r, sizeof(*beam->r) * (beam->np + np_new))) ||
          !(beam->pz = SDDS_Realloc(beam->pz, sizeof(*beam->pz) * (beam->np + np_new))) ||
          !(beam->pr = SDDS_Realloc(beam->pr, sizeof(*beam->pr) * (beam->np + np_new))) ||
          !(beam->r_pphi = SDDS_Realloc(beam->r_pphi, sizeof(*beam->r_pphi) * (beam->np + np_new))) ||
          !(beam->vz = SDDS_Realloc(beam->vz, sizeof(*beam->vz) * (beam->np + np_new))) ||
          !(beam->vr = SDDS_Realloc(beam->vr, sizeof(*beam->vr) * (beam->np + np_new))) ||
          !(beam->vphi = SDDS_Realloc(beam->vphi, sizeof(*beam->vphi) * (beam->np + np_new))) ||
          !(beam->z0 = SDDS_Realloc(beam->z0, sizeof(*beam->z0) * (beam->np + np_new))) ||
          !(beam->t0 = SDDS_Realloc(beam->t0, sizeof(*beam->t0) * (beam->np + np_new))) ||
          !(beam->r0 = SDDS_Realloc(beam->r0, sizeof(*beam->r0) * (beam->np + np_new))) ||
          !(beam->rmin = SDDS_Realloc(beam->rmin, sizeof(*beam->rmin) * (beam->np + np_new))) ||
          !(beam->rmax = SDDS_Realloc(beam->rmax, sizeof(*beam->rmax) * (beam->np + np_new))) ||
          !(beam->zmin = SDDS_Realloc(beam->zmin, sizeof(*beam->zmin) * (beam->np + np_new))) ||
          !(beam->zmax = SDDS_Realloc(beam->zmax, sizeof(*beam->zmax) * (beam->np + np_new))) ||
          !(beam->generation = SDDS_Realloc(beam->generation, sizeof(*beam->generation) * (beam->np + np_new))) ||
          !(beam->particleID = SDDS_Realloc(beam->particleID, sizeof(*beam->particleID) * (beam->np + np_new))) ||
          !(beam->status = SDDS_Realloc(beam->status, sizeof(*beam->status) * (beam->np + np_new))))
        SDDS_Bomb("Memory allocation failure");
      beam->max_np = beam->np + np_new;
#if DEBUG
      fprintf(stderr, "Realloc'd arrays to size %ld from %ld\n", beam->np + np_new, beam->np);
#endif
      if (!(q = SDDS_GetColumn(&SDDSin, "q")) ||
          !(z = SDDS_GetColumn(&SDDSin, "z")) ||
          !(r = SDDS_GetColumn(&SDDSin, "r")) ||
          !(pz = SDDS_GetColumn(&SDDSin, "pz")) ||
          !(pr = SDDS_GetColumn(&SDDSin, "pr")) ||
          !(pphi = SDDS_GetColumn(&SDDSin, "pphi")) ||
          (hasParticleID && !(particleID = SDDS_GetColumn(&SDDSin, "particleID"))))
        SDDS_Bomb("Problem reading data from SDDS file");
#if DEBUG
      fprintf(stderr, "Got column data from SDDS\n");
#endif

      if (hasParticleID)
        {
          for (ip = 0; ip < np_new; ip++)
            if (beam->nextParticleID <= particleID[ip])
              beam->nextParticleID = particleID[ip] + 1;
        }

      i = beam->np;
      for (ip = 0; ip < np_new; ip++, i++)
        {
          beam->Q[i] = q[ip];
          beam->z[i] = beam->z0[ip] = z[ip];
          beam->r[i] = beam->r0[ip] = r[ip];
          beam->rmin[ip] = beam->zmin[ip] = DBL_MAX;
          beam->rmax[ip] = beam->zmax[ip] = -DBL_MAX;
          beam->pz[i] = pz[ip];
          beam->pr[i] = pr[ip];
          beam->t0[i] = -1;
          if (beam->r[i] < 0)
            {
              beam->r[i] *= -1;
              beam->pr[i] *= -1;
            }
          beam->r_pphi[i] = pphi[ip] * r[ip];
          setParticleVelocities(&beam->vz[i], &beam->vr[i], &beam->vphi[i], &gamma,
                                beam->pz[i], beam->pr[i], beam->r_pphi[i], beam->r[i]);
          beam->status[i] = 0;
          beam->generation[i] = 0;
          if (hasParticleID)
            beam->particleID[i] = particleID[ip];
          else
            beam->particleID[i] = beam->nextParticleID++;
          if (beam->z[i] >= EM_problem->zmin && beam->z[i] <= EM_problem->zmax &&
              beam->r[i] < EM_problem->rmax)
            {
              beam->npActive++;
              beam->status[i] = PART_ACTIVE;
            }
        }
#if DEBUG
      fprintf(stderr, "Copied data to arrays\n");
#endif
      free(q);
      free(z);
      free(r);
      free(pz);
      free(pr);
      free(pphi);
      if (hasParticleID)
        free(particleID);
      beam->np = i;
    }
  if (!SDDS_Terminate(&SDDSin))
    SDDS_Bomb("Problem terminating particle file");
}
