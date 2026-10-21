/*************************************************************************\
 * Copyright (c) 2002 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2002 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: integrate.c
 * purpose: do integration
 *
 * Michael Borland, 1992 
 */
#include "spiffe.h"
#include "space_charge.h"
#ifdef _OPENMP
 #include <omp.h>
#endif

static double *JPartial = NULL;
static long JPartialSize = 0;

static double *getJPartialBuffer(long nThreads, long nCells)
{
  long size;
  size = nThreads * nCells;
  if (size > JPartialSize)
    {
      JPartial = trealloc(JPartial, sizeof(*JPartial) * size);
      JPartialSize = size;
    }
  return JPartial;
}

void setup_space_charge(
                        FIELDS *EM_problem,
                        NAMELIST_TEXT *nl_text /* unparsed parameters in namelist text format */
                        )
{

  radial_sharing = 0;
  longitudinal_sharing = 1;

  process_namelist(&space_charge, nl_text);
  print_namelist(stdout, &space_charge);
  bomb("the space_charge namelist is no longer used--use space_charge=1 in the integrate namelist", NULL);

  EM_problem->modes |= FL_SPACE_CHARGE;
}

void add_particle_currents(
                           FIELDS *EM_problem,
                           BEAM *beam)
{
  long nz, nr;
  double dz, dr;
  double zmin, one_over_piDzDr2;
  /*double zmax, rmax;*/
  long nThreads, nCells, nrStride;
  double *JPart;

  nz = EM_problem->nz;
  nr = EM_problem->nr;
  dz = EM_problem->dz;
  dr = EM_problem->dr;
  zmin = EM_problem->zmin;
  /*zmax = EM_problem->zmax;
    rmax = EM_problem->rmax;*/
  one_over_piDzDr2 = 1 / (PI * dz * sqr(dr));

  /* Each thread scatters into a private copy of the grid, which is summed into
   * the real grid afterwards.  The copies have one extra row and column so that
   * the fp==0 smearing writes at iz==nz or ir==nr stay in bounds.
   */
  nrStride = nr + 2;
  nCells = (nz + 2) * nrStride;
  nThreads = spiffeThreads < 1 ? 1 : spiffeThreads;
  JPart = getJPartialBuffer(nThreads, nCells);

  if (EM_problem->modes & FL_TM_FIELDS)
    {
      if (!EM_problem->Jz)
        bomb("Jz array is unexpectedly NULL in add_particle_currents", NULL);
      if (!EM_problem->Jr)
        bomb("Jr array is unexpectedly NULL in add_particle_currents", NULL);
      /* add Jz currents */
#pragma omp parallel num_threads(nThreads)
      {
      long ip, iz, ir, it;
      double fp, fm, Jm, Jp, Q_over_piDzDr2, QVG;
      double *myJ;
#ifdef _OPENMP
      it = omp_get_thread_num();
#else
      it = 0;
#endif
      myJ = JPart + it * nCells;
      memset(myJ, 0, sizeof(*myJ) * nCells);
#pragma omp for schedule(static)
      for (ip = 0; ip < beam->np; ip++)
        {
          if (!(beam->status[ip] & PART_ACTIVE))
            continue;
          if (EM_problem->modes & FL_RADIAL_SMEAR)
            /* compute the ir index for the Jz cell below the particle */
            ir = beam->r[ip] / dr;
          else
            /* compute the ir index for the Jz cell closest to the particle */
            ir = beam->r[ip] / dr + 0.5;
          if (EM_problem->modes & FL_LONGIT_SMEAR)
            /* compute the iz index for the Jz cell to the left of the particle */
            iz = (beam->z[ip] - zmin) / dz + 0.5;
          else
            /* compute the iz index for the Jz cell closest to the particle */
            iz = (beam->z[ip] - zmin) / dz + 1.0;
          if (iz < 0 || iz > nz)
            continue;
          if (ir < 0 || ir > nr)
            continue;
          Q_over_piDzDr2 = beam->Q[ip] * one_over_piDzDr2;
          QVG = beam->vz[ip] * Q_over_piDzDr2;
          if (EM_problem->modes & FL_RADIAL_SMEAR)
            /* determine weight for points at top and bottom of grid square */
            fm = 1 - (fp = beam->r[ip] / dr - ir);
          else
            {
              fm = 1;
              fp = 0;
            }
          /* compute the current densities at top and bottom using volume factors */
          if (ir == 0)
            {
              Jm = 4 * fm * QVG;
              Jp = fp * QVG / (sqr(1.5) - sqr(0.5));
            }
          else if (ir < nr - 1)
            {
              Jm = fm * QVG / (sqr(ir + 0.5) - sqr(ir - 0.5));
              Jp = fp * QVG / (sqr(ir + 1.5) - sqr(ir + 0.5));
            }
          else
            {
              Jm = 0;
              Jp = fp * QVG / (sqr(ir + 0.5) - sqr(ir - 0.5));
              ir--;
            }
          if (EM_problem->modes & FL_LONGIT_SMEAR)
            {
              /* determine weights for grid points to left and right */
              if (iz == 0)
                {
                  fp = 1;
                  fm = 0;
                }
              else if (iz == nz)
                {
                  fm = 1;
                  fp = 0;
                }
              else
                fm = 1 - (fp = (beam->z[ip] - zmin) / dz + 0.5 - iz);
            }
          else
            {
              fm = 1;
              fp = 0;
            }
          /* distribute current density to grid points */
          myJ[iz * nrStride + ir] += fm * Jm;
          if (EM_problem->modes & FL_LONGIT_SMEAR)
            myJ[(iz + 1) * nrStride + ir] += fp * Jm;
          if (EM_problem->modes & FL_RADIAL_SMEAR)
            {
              myJ[iz * nrStride + ir + 1] += fm * Jp;
              if (EM_problem->modes & FL_LONGIT_SMEAR)
                myJ[(iz + 1) * nrStride + ir + 1] += fp * Jp;
            }
        }
      /* implicit barrier above: all scatters are complete */
#pragma omp for schedule(static)
      for (iz = 0; iz <= nz; iz++)
        {
          double *Jrow = EM_problem->Jz[iz];
          for (it = 0; it < nThreads; it++)
            {
              double *prow = JPart + it * nCells + iz * nrStride;
              for (ir = 0; ir <= nr; ir++)
                Jrow[ir] += prow[ir];
            }
        }
      }

      /* add Jr currents */
#pragma omp parallel num_threads(nThreads)
      {
      long ip, iz, ir, it;
      double fp, fm, Jm, Jp, Q_over_piDzDr2, QVG;
      double *myJ;
#ifdef _OPENMP
      it = omp_get_thread_num();
#else
      it = 0;
#endif
      myJ = JPart + it * nCells;
      memset(myJ, 0, sizeof(*myJ) * nCells);
#pragma omp for schedule(static)
      for (ip = 0; ip < beam->np; ip++)
        {
          if (!(beam->status[ip] & PART_ACTIVE))
            continue;
          if (EM_problem->modes & FL_RADIAL_SMEAR)
            /* compute the ir index for the Jr cell below the particle */
            ir = beam->r[ip] / dr + 0.5;
          else
            /* compute the ir index for the Jr cell closest to the particle */
            ir = beam->r[ip] / dr + 1.0;
          if (EM_problem->modes & FL_LONGIT_SMEAR)
            /* compute the iz index for the Jr cell to the left of the particle */
            iz = (beam->z[ip] - zmin) / dz;
          else
            /* compute the iz index for the Jr cell closest to the particle */
            iz = (beam->z[ip] - zmin) / dz + 0.5;
          if (iz < 0 || iz > nz)
            continue;
          if (ir < 0 || ir >= nr)
            continue;
          Q_over_piDzDr2 = beam->Q[ip] * one_over_piDzDr2;
          if (EM_problem->modes & FL_RADIAL_SMEAR)
            /* determine weight for points at top and bottom */
            fm = 1 - (fp = beam->r[ip] / dr + 0.5 - ir);
          else
            {
              fm = 1;
              fp = 0;
            }
          QVG = beam->vr[ip] * Q_over_piDzDr2;
          /* compute the current densities at top and bottom using volume factors */
          if (ir == 0)
            {
              Jm = 0;
              Jp = QVG;
            }
          else if (ir < nr - 1)
            {
              Jm = fm * QVG / (sqr(ir) - sqr(ir - 1));
              Jp = fp * QVG / (sqr(ir + 1) - sqr(ir));
            }
          else
            {
              Jm = 0;
              ir--;
              Jp = QVG / (sqr(ir + 1) - sqr(ir));
            }
          if (EM_problem->modes & FL_LONGIT_SMEAR)
            {
              /* determine weight for points at left and right of grid square */
              if (iz == nz)
                {
                  fm = 1;
                  fp = 0;
                }
              else
                fm = 1 - (fp = (beam->z[ip] - zmin) / dz - iz);
            }
          else
            {
              fm = 1;
              fp = 0;
            }
          /* distribute current density to grid points */
          myJ[iz * nrStride + ir] += fm * Jm;
          if (EM_problem->modes & FL_LONGIT_SMEAR)
            myJ[(iz + 1) * nrStride + ir] += fp * Jm;
          if (EM_problem->modes & FL_RADIAL_SMEAR)
            {
              myJ[iz * nrStride + ir + 1] += fm * Jp;
              if (EM_problem->modes & FL_LONGIT_SMEAR)
                myJ[(iz + 1) * nrStride + ir + 1] += fp * Jp;
            }
        }
      /* implicit barrier above: all scatters are complete */
#pragma omp for schedule(static)
      for (iz = 0; iz <= nz; iz++)
        {
          double *Jrow = EM_problem->Jr[iz];
          for (it = 0; it < nThreads; it++)
            {
              double *prow = JPart + it * nCells + iz * nrStride;
              for (ir = 0; ir <= nr; ir++)
                Jrow[ir] += prow[ir];
            }
        }
      }
    }

  if (EM_problem->modes & FL_TE_FIELDS)
    {
      if (!EM_problem->Jphi)
        bomb("Jphi array is unexpectedly NULL in add_particle_currents", NULL);
      /* add Jphi currents */
#pragma omp parallel num_threads(nThreads)
      {
      long ip, iz, ir, it;
      double fp, fm, Jm, Jp, Q_over_piDzDr2, QVG;
      double *myJ;
#ifdef _OPENMP
      it = omp_get_thread_num();
#else
      it = 0;
#endif
      myJ = JPart + it * nCells;
      memset(myJ, 0, sizeof(*myJ) * nCells);
#pragma omp for schedule(static)
      for (ip = 0; ip < beam->np; ip++)
        {
          if (!(beam->status[ip] & PART_ACTIVE))
            continue;
          if (EM_problem->modes & FL_RADIAL_SMEAR)
            /* compute the ir index for the Jphi cell below the particle */
            ir = beam->r[ip] / dr;
          else
            /* compute the ir index for the Jphi cell closest to the particle */
            ir = beam->r[ip] / dr + 0.5;
          if (EM_problem->modes & FL_LONGIT_SMEAR)
            /* compute the iz index for the Jphi cell to the left of the particle */
            iz = (beam->z[ip] - zmin) / dz;
          else
            /* compute the iz index for the Jphi cell closest to the particle */
            iz = (beam->z[ip] - zmin) / dz + 0.5;
          if (iz < 0 || iz > nz)
            continue;
          if (ir < 0 || ir > nr)
            continue;
          Q_over_piDzDr2 = beam->Q[ip] * one_over_piDzDr2;
          if (EM_problem->modes & FL_RADIAL_SMEAR)
            /* determine weight for points at top and bottom of grid square */
            fm = 1 - (fp = beam->r[ip] / dr - ir);
          else
            {
              fm = 1;
              fp = 0;
            }
          QVG = beam->vphi[ip] * Q_over_piDzDr2;
          /* compute the current densities at top and bottom using volume factors */
          if (ir == 0)
            {
              Jm = 4 * fm * QVG;
              Jp = fp * QVG / (sqr(1.5) - sqr(0.5));
            }
          else if (ir < nr - 1)
            {
              Jm = fm * QVG / (sqr(ir + 0.5) - sqr(ir - 0.5));
              Jp = fp * QVG / (sqr(ir + 1.5) - sqr(ir + 0.5));
            }
          else
            {
              Jm = 0;
              Jp = fp * QVG / (sqr(ir + 0.5) - sqr(ir - 0.5));
              ir--;
            }
          if (EM_problem->modes & FL_LONGIT_SMEAR)
            {
              /* determine weights for grid points to left and right */
              if (iz == nz - 1)
                fm = 1 - (fp = 0);
              else
                fm = 1 - (fp = (beam->z[ip] - zmin) / dz - iz);
            }
          else
            {
              fm = 1;
              fp = 0;
            }
          /* distribute current density to grid points */
          myJ[iz * nrStride + ir] += fm * Jm;
          if (EM_problem->modes & FL_LONGIT_SMEAR)
            myJ[(iz + 1) * nrStride + ir] += fp * Jm;
          if (EM_problem->modes & FL_RADIAL_SMEAR)
            {
              myJ[iz * nrStride + ir + 1] += fm * Jp;
              if (EM_problem->modes & FL_LONGIT_SMEAR)
                myJ[(iz + 1) * nrStride + ir + 1] += fp * Jp;
            }
        }
      /* implicit barrier above: all scatters are complete */
#pragma omp for schedule(static)
      for (iz = 0; iz <= nz; iz++)
        {
          double *Jrow = EM_problem->Jphi[iz];
          for (it = 0; it < nThreads; it++)
            {
              double *prow = JPart + it * nCells + iz * nrStride;
              for (ir = 0; ir <= nr; ir++)
                Jrow[ir] += prow[ir];
            }
        }
      }
    }
}
