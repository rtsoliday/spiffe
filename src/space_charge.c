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
  long ip;
  long iz, ir, nz, nr;
  double fp, fm, Jm, Jp, dz, dr;
  double Q_over_piDzDr2, zmin, one_over_piDzDr2;
  /*double zmax, rmax;*/
  double QVG;

  nz = EM_problem->nz;
  nr = EM_problem->nr;
  dz = EM_problem->dz;
  dr = EM_problem->dr;
  zmin = EM_problem->zmin;
  /*zmax = EM_problem->zmax;
    rmax = EM_problem->rmax;*/
  one_over_piDzDr2 = 1 / (PI * dz * sqr(dr));

  if (EM_problem->modes & FL_TM_FIELDS)
    {
      if (!EM_problem->Jz)
        bomb("Jz array is unexpectedly NULL in add_particle_currents", NULL);
      if (!EM_problem->Jr)
        bomb("Jr array is unexpectedly NULL in add_particle_currents", NULL);
      /* add Jz currents */
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
          EM_problem->Jz[iz][ir] += fm * Jm;
          if (EM_problem->modes & FL_LONGIT_SMEAR)
            EM_problem->Jz[iz + 1][ir] += fp * Jm;
          if (EM_problem->modes & FL_RADIAL_SMEAR)
            {
              EM_problem->Jz[iz][ir + 1] += fm * Jp;
              if (EM_problem->modes & FL_LONGIT_SMEAR)
                EM_problem->Jz[iz + 1][ir + 1] += fp * Jp;
            }
        }

      /* add Jr currents */
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
          EM_problem->Jr[iz][ir] += fm * Jm;
          if (EM_problem->modes & FL_LONGIT_SMEAR)
            EM_problem->Jr[iz + 1][ir] += fp * Jm;
          if (EM_problem->modes & FL_RADIAL_SMEAR)
            {
              EM_problem->Jr[iz][ir + 1] += fm * Jp;
              if (EM_problem->modes & FL_LONGIT_SMEAR)
                EM_problem->Jr[iz + 1][ir + 1] += fp * Jp;
            }
        }
    }

  if (EM_problem->modes & FL_TE_FIELDS)
    {
      if (!EM_problem->Jphi)
        bomb("Jphi array is unexpectedly NULL in add_particle_currents", NULL);
      /* add Jphi currents */
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
          EM_problem->Jphi[iz][ir] += fm * Jm;
          if (EM_problem->modes & FL_LONGIT_SMEAR)
            EM_problem->Jphi[iz + 1][ir] += fp * Jm;
          if (EM_problem->modes & FL_RADIAL_SMEAR)
            {
              EM_problem->Jphi[iz][ir + 1] += fm * Jp;
              if (EM_problem->modes & FL_LONGIT_SMEAR)
                EM_problem->Jphi[iz + 1][ir + 1] += fp * Jp;
            }
        }
    }
}
