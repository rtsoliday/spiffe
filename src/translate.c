/*************************************************************************\
 * Copyright (c) 2002 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2002 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: translate.c
 * purpose: do translation of problem in z
 *
 * Michael Borland, 1992 
 */
#include "spiffe.h"
#include "translate.h"

static long do_translation = 0;
static long region_reduced;

void setup_translation(
                       FIELDS *EM_problem,
                       NAMELIST_TEXT *nl_text /* unparsed parameters in namelist text format */
                       )
{

  z_trigger = EM_problem->zmax - 5 * EM_problem->dz;
  z_lower = EM_problem->zmin;
  z_upper = EM_problem->zmax;

  process_namelist(&translate, nl_text);
  print_namelist(stdout, &translate);

  if (z_trigger < EM_problem->zmin + EM_problem->dz || z_trigger > EM_problem->zmax - EM_problem->dz)
    bomb("z_trigger is not in range [zmin+dz, zmax-dz]", NULL);
  do_translation = 1;
  region_reduced = 0;
}

void translate_problem(BEAM *beam, FIELDS *EM_problem, ANTENNA *antenna, long n_antennae, RESISTOR *resistor)
{
  long iz, ir, ip, i;
  long nz, nr, diz;
  double *ptr[10];

  if (!do_translation)
    return;

  for (ip = 0; ip < beam->np; ip++)
    if (beam->z[ip] > z_trigger)
      break;
  if (ip == beam->np)
    {
      return;
    }

  if (!region_reduced && z_lower != z_upper && !(z_lower <= EM_problem->zmin && z_upper >= EM_problem->zmax))
    {
      /* reduce the size of the region in which fields are integrated */
      /* should probably free memory here instead of just writing over the pointers */
      nz = EM_problem->nz;
      EM_problem->upper_bc = EM_problem->right_bc =
        EM_problem->left_bc = NEUMANN;
      EM_problem->lower_bc = DIRICHLET;
      if ((diz = (z_lower - EM_problem->zmin) / EM_problem->dz) > 0)
        {
          for (iz = 0; iz <= nz - diz; iz++)
            EM_problem->metal_flags[iz] = EM_problem->metal_flags[iz + diz];
          if (EM_problem->modes & FL_TM_FIELDS)
            for (iz = 0; iz <= nz - diz; iz++)
              {
                EM_problem->Ez[iz] = EM_problem->Ez[iz + diz];
                EM_problem->Jz[iz] = EM_problem->Jz[iz + diz];
                EM_problem->Er[iz] = EM_problem->Er[iz + diz];
                EM_problem->Jr[iz] = EM_problem->Jr[iz + diz];
                EM_problem->Bphi[iz] = EM_problem->Bphi[iz + diz];
                EM_problem->Bphi2[iz] = EM_problem->Bphi2[iz + diz];
              }
          if (EM_problem->modes & FL_TE_FIELDS)
            for (iz = 0; iz <= nz - diz; iz++)
              {
                EM_problem->Ephi[iz] = EM_problem->Ephi[iz + diz];
                EM_problem->Jphi[iz] = EM_problem->Jphi[iz + diz];
                EM_problem->Br[iz] = EM_problem->Br[iz + diz];
                EM_problem->Br2[iz] = EM_problem->Br2[iz + diz];
                EM_problem->Bz[iz] = EM_problem->Bz[iz + diz];
                EM_problem->Bz2[iz] = EM_problem->Bz2[iz + diz];
              }
          EM_problem->nz -= diz;
          for (i = 0; i < n_antennae; i++)
            {
              if (antenna[i].direction == 0)
                {
                  antenna[i].istart -= diz;
                  antenna[i].iend -= diz;
                }
              else
                antenna[i].iposition -= diz;
            }
          for (i = 0; i < resistor->n_resistors; i++)
            {
              if (resistor->direction[i] == 0)
                {
                  resistor->istart[i] -= diz;
                  resistor->iend[i] -= diz;
                }
              else
                resistor->iposition[i] -= diz;
            }
          EM_problem->zmin += diz * EM_problem->dz;
        }
      if ((diz = (EM_problem->zmax - z_upper)) > 0)
        {
          EM_problem->nz -= diz;
          EM_problem->zmax -= diz * EM_problem->dz;
        }
      printf("\nField integration region reduced to z:[%e, %e] (%" PRId32 " mesh lines)\n",
             EM_problem->zmin, EM_problem->zmax, EM_problem->nz);
      if (fabs((EM_problem->zmax - EM_problem->zmin) / (EM_problem->nz - 1) - EM_problem->dz) > 1e-6)
        bomb("region reduction changed dz unacceptably!", NULL);
      region_reduced = 1;
    }

  /* translate the field pointers */
  nz = EM_problem->nz;
  nr = EM_problem->nr;
  if (EM_problem->modes & FL_TM_FIELDS)
    {
      ptr[0] = EM_problem->Ez[0];
      ptr[1] = EM_problem->Jz[0];
      ptr[2] = EM_problem->Er[0];
      ptr[3] = EM_problem->Jr[0];
      ptr[4] = EM_problem->Bphi[0];
      ptr[5] = EM_problem->Bphi2[0];
      for (iz = 0; iz < nz; iz++)
        {
          EM_problem->Ez[iz] = EM_problem->Ez[iz + 1];
          EM_problem->Jz[iz] = EM_problem->Jz[iz + 1];
          EM_problem->Er[iz] = EM_problem->Er[iz + 1];
          EM_problem->Jr[iz] = EM_problem->Jr[iz + 1];
          EM_problem->Bphi[iz] = EM_problem->Bphi[iz + 1];
          EM_problem->Bphi2[iz] = EM_problem->Bphi2[iz + 1];
        }
      EM_problem->Ez[nz] = ptr[0];
      EM_problem->Jz[nz] = ptr[1];
      EM_problem->Er[nz] = ptr[2];
      EM_problem->Jr[nz] = ptr[3];
      EM_problem->Bphi[nz] = ptr[4];
      EM_problem->Bphi2[nz] = ptr[5];
      for (ir = 0; ir <= nr; ir++)
        EM_problem->Ez[nz][ir] =
          EM_problem->Jz[nz][ir] =
          EM_problem->Er[nz][ir] =
          EM_problem->Jr[nz][ir] =
          EM_problem->Bphi[nz][ir] =
          EM_problem->Bphi2[nz][ir] = 0;
    }
  if (EM_problem->modes & FL_TE_FIELDS)
    {
      ptr[0] = EM_problem->Ephi[0];
      ptr[1] = EM_problem->Jphi[0];
      ptr[2] = EM_problem->Br[0];
      ptr[3] = EM_problem->Br2[0];
      ptr[4] = EM_problem->Bz[0];
      ptr[5] = EM_problem->Bz2[0];
      for (iz = 0; iz < nz; iz++)
        {
          EM_problem->Ephi[iz] = EM_problem->Ephi[iz + 1];
          EM_problem->Jphi[iz] = EM_problem->Jphi[iz + 1];
          EM_problem->Br[iz] = EM_problem->Br[iz + 1];
          EM_problem->Br2[iz] = EM_problem->Br2[iz + 1];
          EM_problem->Bz[iz] = EM_problem->Bz[iz + 1];
          EM_problem->Bz2[iz] = EM_problem->Bz2[iz + 1];
        }
      EM_problem->Ephi[nz] = ptr[0];
      EM_problem->Jphi[nz] = ptr[1];
      EM_problem->Br[nz] = ptr[2];
      EM_problem->Br2[nz] = ptr[3];
      EM_problem->Bz[nz] = ptr[4];
      EM_problem->Bz2[nz] = ptr[5];
      for (ir = 0; ir <= nr; ir++)
        EM_problem->Ephi[nz][ir] =
          EM_problem->Jphi[nz][ir] =
          EM_problem->Br[nz][ir] =
          EM_problem->Br2[nz][ir] =
          EM_problem->Bz[nz][ir] =
          EM_problem->Bz2[nz][ir] = 0;
    }

  EM_problem->zmin += EM_problem->dz;
  EM_problem->zmax += EM_problem->dz;
  EM_problem->zt += EM_problem->dz;
  z_trigger += EM_problem->dz;

  printf("problem translated by %em for total translation of %em at time %es\n",
         EM_problem->dz, EM_problem->zt, EM_problem->time);
}
