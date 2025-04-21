/*************************************************************************\
 * Copyright (c) 2002 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2002 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: advance_particles.
 * purpose: advance positions of electrons by one time step
 *
 * Michael Borland, 1992
 */
#include "spiffe.h"

/* for non-ANSI systems */
#if defined(vxWorks) || defined(__BORLANDC__) || defined(linux) || defined(_AIX) || defined(__rtems__) || defined(_MINGW)
#  ifndef HUGE_VAL
#    define HUGE_VAL HUGE
#  else
#    ifndef HUGE
#      define HUGE HUGE_VAL
#    endif
#  endif
#endif

short checkLossCode(long iz, long ir, FIELDS *EM_problem)
{
  long lossCode, iz1, ir1;
  lossCode = 0;
  if (EM_problem->metal_flags[iz][ir] & FL_IS_INSIDE)
    {
      lossCode = 3;
    }
  else if (EM_problem->metal_flags[iz][ir] & FL_IS_SURFACE)
    {
      iz1 = iz + 1;
      if (iz1 > EM_problem->nz)
        {
          lossCode = 4;
        }
      else
        {
          ir1 = ir + 1;
          if (ir1 > EM_problem->nr)
            {
              lossCode = 5;
            }
          else
            {
              if (EM_problem->metal_flags[iz1][ir1] & FL_IS_METAL &&
                  EM_problem->metal_flags[iz][ir1] & FL_IS_METAL &&
                  EM_problem->metal_flags[iz1][ir] & FL_IS_METAL)
                {
                  lossCode = 6;
                }
            }
        }
    }
  return lossCode;
}

void advance_particles(
                       BEAM *beam,
                       FIELDS *EM_problem,
                       TEST_BEAM_DEFN *test_beam_defn,
                       SCREEN *screen_defn,
                       long n_screens,
                       SDDS_DATASET *SDDSlost,
                       double zForcesStart,
                       SECONDARY_EMISSION *secondaryEmission,
                       long nSecondaryDefns)
{
  /*double fa, fb, fr, fl;*/
  double z0, r0, z, r, pz = 0, pr = 0, vz = 0, vr = 0, r_pphi = 0, vphi = 0, gamma, dtp, pphi;
  double pz0 = 0, pr0 = 0, r_pphi0 = 0, pphi0 = 0, r0Save, z0Save;
  long ip, is, ib, iz, ir, r_was_negative, is_lost, i;
  //long n_lost;
  /*long longit_interp, radial_interp, nr, nz;*/
  short lossCode;
  double zmin, dz, dr, zmax, rmax;
  double Ez, Er, Bphi, Ephi, Br, Bz, dt;
  /*double EzOOE, ErOOE, BphiOOE;*/
  double z_tb1 = 0, r_tb1 = 0, adjustFactor;
  /*double dtAdjust;*/
  long maxLost, nExtra = 0, newLost = 0;
  double max_r, ave_z, ave_vz;
  short directionCode;
#if defined(DEBUG)
  long izNew, irNew;
#endif

  if (N_SCREEN_QUANS != 8)
    bomb("problem with programming of screens", NULL);

  dt = EM_problem->dti;
  dz = EM_problem->dz;
  dr = EM_problem->dr;
  zmin = EM_problem->zmin;
  zmax = EM_problem->zmax;
  rmax = EM_problem->rmax;
  /*nz = EM_problem->nz;
    nr = EM_problem->nr;*/
  //n_lost = 0;
  maxLost = 0;
  /*longit_interp = EM_problem->modes&FL_LONGIT_INTERP;
    radial_interp = EM_problem->modes&FL_RADIAL_INTERP;*/

  if (test_beam_defn->electrons_per_macroparticle)
    {
      r_tb1 = test_beam_defn->r_start_deceleration;
      z_tb1 = test_beam_defn->z_start_deceleration;
    }

  for (is = 0; is < n_screens; is++)
    screen_defn[is].nb = 0;
  max_r = 0;
  ave_z = ave_vz = 0;
  if (beam->np == 0)
    nExtra = 0;
  EM_problem->EzOOEMin = HUGE;
  EM_problem->EzOOEMax = -HUGE;

#if DEBUG
  fprintf(stdout, "%ld particles, %ld active particles, %ld extra particles\n",
          beam->np, beam->npActive, nExtra);
  fflush(stdout);
#endif
  for (ip = 0; ip < beam->np + nExtra; ip++)
    {
#if DEBUG
      fprintf(stdout, "Particle %ld of %ld: %s, z=%e, r=%e, pz=%e, pr=%e, r_pphi=%e, vz=%e, vr=%e, vphi=%e\n",
              ip, beam->np+nExtra, beam->status[ip]&PART_ACTIVE?"active":"inactive",
              beam->z[ip], beam->r[ip], beam->pz[ip], beam->pr[ip], beam->r_pphi[ip],
              beam->vz[ip], beam->vr[ip], beam->vphi[ip]);
      fflush(stdout);
#endif

      /*dtAdjust = 0;*/
      adjustFactor = 1;
      z = z0 = r = r0 = z0Save = r0Save = 0;

      if (ip < beam->np)
        {

          /* (z, r) values are at t-dt (synchronized with Ez, Er). */
          z = z0 = z0Save = beam->z[ip];
          r = r0 = r0Save = beam->r[ip];

          /* update statistics on particle motion */
          if (z > beam->zmax[ip])
            beam->zmax[ip] = z;
          if (z < beam->zmin[ip])
            beam->zmin[ip] = z;
          if (r > beam->rmax[ip])
            beam->rmax[ip] = r;
          if (r < beam->rmin[ip])
            beam->rmin[ip] = r;

          /* (vz, vr, pz, pr) values are at t-dt/2 (synchronized with Bphi). */
          vz = beam->vz[ip];
          vr = beam->vr[ip];

          /* use (vz, vr) to advance (z, r) to time t */
          z += vz * dt;
          r += vr * dt;
          if (!(beam->status[ip] & PART_ACTIVE))
            {
              iz = (z - zmin) / dz;
              ir = fabs(r) / dr;
              if (!(!(z >= EM_problem->zmin && z <= EM_problem->zmax && r <= EM_problem->rmax) ||
                    checkLossCode(iz, ir, EM_problem))) {
                beam->status[ip] |= PART_ACTIVE;
                beam->npActive++;
                beam->t0[ip] = EM_problem->time;
                beam->r0[ip] = beam->r[ip];
                beam->z0[ip] = beam->z[ip];
                /* Assume the particle entered on the left. This factor is used to adjust the integral 
                 * of the force applied on this step 
                 */
                if (iz > 0 && checkLossCode(iz - 1, ir, EM_problem))
                  adjustFactor = (z - (iz * EM_problem->dz + EM_problem->zmin)) / (z - z0);
                else
                  adjustFactor = (z - EM_problem->zmin) / (z - z0);
                if (adjustFactor > 1)
                  adjustFactor = 1;
              }
            }
          if (zForcesStart > EM_problem->zmin && zForcesStart < EM_problem->zmax)
            {
              if (z0 < zForcesStart)
                {
                  if (z < zForcesStart)
                    adjustFactor = 0;
                  else
                    {
                      adjustFactor = (z - zForcesStart) / (z - z0);
                    }
                }
            }
          pz = pz0 = beam->pz[ip];
          pr = pr0 = beam->pr[ip];
          r_pphi = r_pphi0 = beam->r_pphi[ip];
          vphi = beam->vphi[ip];
          if (r > max_r)
            max_r = r;
          ave_z += z;
          ave_vz += vz;
        }
      else
        {
          /* set values for field interpolation test output */
          if (ip == beam->np)
            {
              ave_z /= beam->npActive;
              ave_vz /= beam->npActive;
            }
          z = z0 = ave_z;
          r = r0 = (ip - beam->np) * (max_r / nExtra);
          vz = ave_vz;
          pz = pz0 = vz / c_mks / sqrt(1 - sqr(vz / c_mks));
          pr = pr0 = 0;
          vr = vphi = 0;
          r_pphi = r_pphi0 = 0;
        }

      is_lost = 0;
      lossCode = 0;
      if ((beam->status[ip] & PART_ACTIVE) && (z < zmin || z > zmax || r > rmax))
        {
          lossCode = 1;
          is_lost = 1;
        }

      /* to determine if particle is inside metal, compute indices of left/bottom grid point */
      iz = (z - zmin) / dz;
      ir = fabs(r) / dr;
#if DEBUG
      izNew = iz;
      irNew = ir;
#endif
      if ((beam->status[ip] & PART_ACTIVE) && (iz > EM_problem->nz || ir > EM_problem->nr || iz < 0 || ir < 0))
        {
          lossCode = 2;
          is_lost = 1;
        }

      if (!is_lost && (beam->status[ip] & PART_ACTIVE))
        {
          if ((lossCode = checkLossCode(iz, ir, EM_problem)))
            is_lost = 1;
        }

#if DEBUG
      if (is_lost)
        fprintf(stdout, "particle %ld lost: z=%e, r=%e\n", ip, z, r);
#endif
      r_was_negative = 0;
      if (r < 0)
        {
          r *= -1;
          pr *= -1;
          r_was_negative = 1;
        }

      if (!is_lost && (beam->status[ip] & PART_ACTIVE))
        {
          /* update statistics on particle motion */
          if (z > beam->zmax[ip])
            beam->zmax[ip] = z;
          if (z < beam->zmin[ip])
            beam->zmin[ip] = z;
          if (r > beam->rmax[ip])
            beam->rmax[ip] = r;
          if (r < beam->rmin[ip])
            beam->rmin[ip] = r;
        }

      /* compute indices of closest grid point */
      iz = (z - zmin) / dz + 0.5;
      ir = r / dr + 0.5;

      if (!(test_beam_defn->electrons_per_macroparticle &&
            (test_beam_defn->z_force || test_beam_defn->r_force)))
        {
          if (beam->status[ip] & PART_ACTIVE) {
            interpolateFields(&Ez, &Er, &Ephi, &Bz, &Br, &Bphi, z, r, EM_problem);
          } else {
            /* calculate fields at new position */
            Ez = Er = Ephi = Bz = Br = Bphi = 0;
          }

          if (test_beam_defn->electrons_per_macroparticle)
            beam->stiffness = 1;
          /* advance (pr, pz, r_pphi) from t-dt/2 to t+dt/2 */
          pphi = (r ? r_pphi / r : 0);
          gamma = sqrt(sqr(pr) + sqr(pz) + sqr(pphi) + 1);
          pr += adjustFactor * dt * ((r ? sqr(pphi) * c_mks / gamma / r : 0) + beam->QoMC * (Er + vphi * Bz - vz * Bphi)) / beam->stiffness;
          pz += adjustFactor * dt * beam->QoMC * (Ez + vr * Bphi - vphi * Br) / beam->stiffness;
          r_pphi += adjustFactor * dt * beam->QoMC * r * (Ephi + vz * Br - vr * Bz) / beam->stiffness;

        }
      else
        {
          if (ip >= beam->np)
            continue;
          /* test particles with predefined acceleration: no effect from fields */
          pphi = (r ? r_pphi / r : 0);
          gamma = sqrt(sqr(pr) + sqr(pz) + sqr(pphi) + 1);
          pr += dt * ((r ? sqr(pphi) * c_mks / gamma / r : 0));
          if (z <= z_tb1)
            pz += dt * test_beam_defn->z_force;
          else
            {
              pz -= dt * test_beam_defn->z_force;
              if (pz < 0)
                pz = 0;
            }
          if (r_was_negative)
            pr = 0;
          else if (r <= r_tb1)
            pr += dt * test_beam_defn->r_force;
          else
            pr -= dt * test_beam_defn->r_force;
        }

      beam->z[ip] = z;
      beam->r[ip] = r;
      beam->pz[ip] = pz;
      beam->pr[ip] = pr;
      beam->r_pphi[ip] = r_pphi;
#if DEBUG
      fprintf(stdout, "Particle %ld coords: z=%e, r=%e, pz=%e, pr=%e, r_pphi=%e\n",
              ip, z, r, pz, pr, r_pphi);
      fflush(stdout);
#endif
      setParticleVelocities(&beam->vz[ip], &beam->vr[ip], &beam->vphi[ip], &gamma,
                            pz, pr, r_pphi, r);
#if DEBUG
      fprintf(stdout, "Velocities: (%e, %e, %e)\n",
              beam->vz[ip], beam->vr[ip], beam->vphi[ip]);
      fflush(stdout);
#endif

      for (is = 0; is < n_screens; is++)
        {
          if (screen_defn[is].z_position <= z && screen_defn[is].z_position > z0 && screen_defn[is].direction == 1)
            {
              if (is_lost)
                {
                  fprintf(stderr, "*** WARNING: a particle is recorded on screen %s that is lost in the same time-step.\n",
                          screen_defn[is].filename);
                  fprintf(stderr, "The data for this particle is unreliable--consider moving the screen to z-dz.\n");
                }
              if (screen_defn[is].nb >= screen_defn[is].max_nb)
                {
                  screen_defn[is].max_nb += 10;
                  for (i = 0; i < N_SCREEN_QUANS; i++)
                    screen_defn[is].buffer[i] =
                      trealloc(screen_defn[is].buffer[i], sizeof(*screen_defn[is].buffer[i]) * screen_defn[is].max_nb);
                  screen_defn[is].generation =
                    trealloc(screen_defn[is].generation, sizeof(*screen_defn[is].generation) * screen_defn[is].max_nb);
                  screen_defn[is].particleID =
                    trealloc(screen_defn[is].particleID, sizeof(*screen_defn[is].particleID) * screen_defn[is].max_nb);
                  screen_defn[is].rmin =
                    trealloc(screen_defn[is].rmin, sizeof(*screen_defn[is].rmin) * screen_defn[is].max_nb);
                  screen_defn[is].rmax =
                    trealloc(screen_defn[is].rmax, sizeof(*screen_defn[is].rmax) * screen_defn[is].max_nb);
                  screen_defn[is].zmin =
                    trealloc(screen_defn[is].zmin, sizeof(*screen_defn[is].zmin) * screen_defn[is].max_nb);
                  screen_defn[is].zmax =
                    trealloc(screen_defn[is].zmax, sizeof(*screen_defn[is].zmax) * screen_defn[is].max_nb);
                }
              dtp = (screen_defn[is].z_position - z0) / vz;
              screen_defn[is].buffer[0][screen_defn[is].nb] = r0 + vr * dtp;
              screen_defn[is].buffer[1][screen_defn[is].nb] = pz0 + (pz - pz0) / dt * dtp;
              screen_defn[is].buffer[2][screen_defn[is].nb] = pr0 + (pr - pr0) / dt * dtp;
              pphi = r ? r_pphi / r : 0;
              pphi0 = r0 ? r_pphi0 / r0 : 0;
              screen_defn[is].buffer[3][screen_defn[is].nb] = pphi0 + (pphi - pphi0) / dt * dtp;
              screen_defn[is].buffer[4][screen_defn[is].nb] = (EM_problem->time - dt + dtp);
              screen_defn[is].buffer[5][screen_defn[is].nb] = beam->Q[ip];
              screen_defn[is].buffer[6][screen_defn[is].nb] = beam->r0[ip];
              screen_defn[is].buffer[7][screen_defn[is].nb] = beam->t0[ip];
              screen_defn[is].rmin[screen_defn[is].nb] = beam->rmin[ip];
              screen_defn[is].rmax[screen_defn[is].nb] = beam->rmax[ip];
              screen_defn[is].zmin[screen_defn[is].nb] = beam->zmin[ip];
              screen_defn[is].zmax[screen_defn[is].nb] = beam->zmax[ip];
              screen_defn[is].generation[screen_defn[is].nb] = beam->generation[ip];
              screen_defn[is].particleID[screen_defn[is].nb] = beam->particleID[ip];
              screen_defn[is].nb += 1;
            }
          else if (screen_defn[is].z_position >= z && screen_defn[is].z_position < z0 && screen_defn[is].direction == -1)
            {
              if (is_lost)
                {
                  fprintf(stderr, "*** WARNING: a particle is recorded on screen %s that is lost in the same time-step.\n",
                          screen_defn[is].filename);
                  fprintf(stderr, "The data for this particle is unreliable--consider moving the screen to z+dz.\n");
                }
              if (screen_defn[is].nb >= screen_defn[is].max_nb)
                {
                  screen_defn[is].max_nb += 10;
                  for (i = 0; i < N_SCREEN_QUANS; i++)
                    screen_defn[is].buffer[i] =
                      trealloc(screen_defn[is].buffer[i], sizeof(*screen_defn[is].buffer[i]) * screen_defn[is].max_nb);
                  screen_defn[is].generation =
                    trealloc(screen_defn[is].generation, sizeof(*screen_defn[is].generation) * screen_defn[is].max_nb);
                  screen_defn[is].particleID =
                    trealloc(screen_defn[is].particleID, sizeof(*screen_defn[is].particleID) * screen_defn[is].max_nb);
                  screen_defn[is].rmin =
                    trealloc(screen_defn[is].rmin, sizeof(*screen_defn[is].rmin) * screen_defn[is].max_nb);
                  screen_defn[is].rmax =
                    trealloc(screen_defn[is].rmax, sizeof(*screen_defn[is].rmax) * screen_defn[is].max_nb);
                  screen_defn[is].zmin =
                    trealloc(screen_defn[is].zmin, sizeof(*screen_defn[is].zmin) * screen_defn[is].max_nb);
                  screen_defn[is].zmax =
                    trealloc(screen_defn[is].zmax, sizeof(*screen_defn[is].zmax) * screen_defn[is].max_nb);
                }
              dtp = (screen_defn[is].z_position - z0) / vz;
              r0 += vr * dtp;
              screen_defn[is].buffer[0][screen_defn[is].nb] = r0;
              screen_defn[is].buffer[1][screen_defn[is].nb] = pz;
              screen_defn[is].buffer[2][screen_defn[is].nb] = pr;
              screen_defn[is].buffer[3][screen_defn[is].nb] = r0 ? r_pphi / r0 : 0;
              screen_defn[is].buffer[4][screen_defn[is].nb] = (EM_problem->time - dt + dtp);
              screen_defn[is].buffer[5][screen_defn[is].nb] = beam->Q[ip];
              screen_defn[is].buffer[6][screen_defn[is].nb] = beam->r0[ip];
              screen_defn[is].buffer[7][screen_defn[is].nb] = beam->t0[ip];
              screen_defn[is].generation[screen_defn[is].nb] = beam->generation[ip];
              screen_defn[is].particleID[screen_defn[is].nb] = beam->particleID[ip];
              screen_defn[is].rmin[screen_defn[is].nb] = beam->rmin[ip];
              screen_defn[is].rmax[screen_defn[is].nb] = beam->rmax[ip];
              screen_defn[is].zmin[screen_defn[is].nb] = beam->zmin[ip];
              screen_defn[is].zmax[screen_defn[is].nb] = beam->zmax[ip];
              screen_defn[is].nb += 1;
            }
        }

#if DEBUG
      fprintf(stdout, "Done processing screens\n");
      fflush(stdout);
#endif

      if (is_lost)
        {
          double p;
          double zLost, rLost;
          p = sqr(beam->pz[ip]) + sqr(beam->pr[ip]);
          if (beam->r[ip])
            p += sqr(beam->r_pphi[ip] / beam->r[ip]);
          p = sqrt(p);
          findLossCoordinates(&zLost, &rLost, &directionCode, z0Save, r0Save, z, r * (r_was_negative ? -1 : 1), EM_problem);
          if (nSecondaryDefns)
            queueSecondaryEmission(secondaryEmission, nSecondaryDefns, p,
                                   zLost, rLost, directionCode, beam->Q[ip],
                                   beam->generation[ip], EM_problem, beam);

          beam->status[ip] &= ~PART_ACTIVE;
          beam->status[ip] |= PART_LOST;
          if (SDDSlost)
            {
              if (newLost && newLost >= maxLost)
                {
                  if (!SDDS_WritePage(SDDSlost))
                    {
                      SDDS_SetError("Problem writing lost particle data.");
                      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                    }
                  newLost = maxLost = 0;
                }
              if (maxLost == 0)
                {
                  if (!SDDS_StartPage(SDDSlost, maxLost = 1000))
                    {
                      SDDS_SetError("Problem writing lost particle data.");
                      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                    }
                }
              if (!SDDS_SetRowValues(SDDSlost, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                                     newLost,
                                     "t", EM_problem->time,
                                     "rLost", (float)rLost, "zLost", (float)zLost,
                                     "pz", (float)beam->pz[ip], "pr", (float)beam->pr[ip],
                                     "r0", (float)beam->r0[ip], "z0", (float)beam->z0[ip],
                                     "q", (float)beam->Q[ip],
                                     "generation", beam->generation[ip],
                                     "particleID", beam->particleID[ip],
#ifdef DEBUG
                                     "lossCode", lossCode,
                                     "iz", izNew, "ir", irNew,
                                     "zLast", z, "rLast", r,
#endif
                                     NULL))
                {
                  SDDS_SetError("Problem writing lost particle data.");
                  SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                }
            }
          beam->np -= 1;
          beam->npActive -= 1;
          beam->Q[ip] = beam->Q[beam->np];
          beam->z[ip] = beam->z[beam->np];
          beam->r[ip] = beam->r[beam->np];
          beam->pz[ip] = beam->pz[beam->np];
          beam->pr[ip] = beam->pr[beam->np];
          beam->vz[ip] = beam->vz[beam->np];
          beam->vr[ip] = beam->vr[beam->np];
          beam->r_pphi[ip] = beam->r_pphi[beam->np];
          beam->vphi[ip] = beam->vphi[beam->np];
          beam->t0[ip] = beam->t0[beam->np];
          beam->r0[ip] = beam->r0[beam->np];
          beam->z0[ip] = beam->z0[beam->np];
          beam->rmin[ip] = beam->rmin[beam->np];
          beam->zmin[ip] = beam->zmin[beam->np];
          beam->rmax[ip] = beam->rmax[beam->np];
          beam->zmax[ip] = beam->zmax[beam->np];
          beam->generation[ip] = beam->generation[beam->np];
          beam->particleID[ip] = beam->particleID[beam->np];
          beam->status[ip] = beam->status[beam->np];
          ip--;
          //n_lost++;
          newLost++;
          continue;
        }
#if DEBUG
      fprintf(stdout, "Done processing lost particles\n");
      fflush(stdout);
#endif
    }

  if (newLost)
    {
      if (SDDSlost && !SDDS_WritePage(SDDSlost))
        {
          SDDS_SetError("Problem writing lost particle data.");
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
        }
      /* printf("%ld particles lost\n", n_lost);  */
    }

  for (is = 0; is < n_screens; is++)
    {
      if (screen_defn[is].nb)
        {
          if (!SDDS_StartTable(&screen_defn[is].table, screen_defn[is].nb))
            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
          for (ib = 0; ib < N_SCREEN_QUANS; ib++)
            if (!SDDS_SetColumn(&screen_defn[is].table, SDDS_SET_BY_INDEX, screen_defn[is].buffer[ib],
                                screen_defn[is].nb, ib))
              SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
          if (!SDDS_SetColumn(&screen_defn[is].table, SDDS_SET_BY_NAME, screen_defn[is].generation,
                              screen_defn[is].nb, "generation") ||
              !SDDS_SetColumn(&screen_defn[is].table, SDDS_SET_BY_NAME, screen_defn[is].particleID,
                              screen_defn[is].nb, "particleID") ||
              !SDDS_SetColumn(&screen_defn[is].table, SDDS_SET_BY_NAME, screen_defn[is].rmin,
                              screen_defn[is].nb, "rMin") ||
              !SDDS_SetColumn(&screen_defn[is].table, SDDS_SET_BY_NAME, screen_defn[is].rmax,
                              screen_defn[is].nb, "rMax") ||
              !SDDS_SetColumn(&screen_defn[is].table, SDDS_SET_BY_NAME, screen_defn[is].zmin,
                              screen_defn[is].nb, "zMin") ||
              !SDDS_SetColumn(&screen_defn[is].table, SDDS_SET_BY_NAME, screen_defn[is].zmax,
                              screen_defn[is].nb, "zMax") ||
              !SDDS_SetParameters(&screen_defn[is].table, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                                  "Time", EM_problem->time, "Timestep", EM_problem->time_step, NULL))
            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
          if (!SDDS_WriteTable(&screen_defn[is].table))
            {
              fprintf(stderr, "Problem writing to screen %s\n", screen_defn[is].filename);
              SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
            }
        }
    }
#if DEBUG
  fprintf(stdout, "Done advancing particles\n");
  fflush(stdout);
#endif
}

void setUpLostParticleFile(SDDS_DATASET *SDDSout, char *filename0, char *rootname)
{
  char *filename;
  filename = compose_filename(filename0, rootname);
  if (!SDDS_InitializeOutput(SDDSout, SDDS_BINARY, 0, NULL, NULL, filename) ||
      !SDDS_DefineSimpleColumn(SDDSout, "t", "s", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "rLost", "m", SDDS_FLOAT) ||
      !SDDS_DefineSimpleColumn(SDDSout, "zLost", "m", SDDS_FLOAT) ||
      !SDDS_DefineSimpleColumn(SDDSout, "pz", NULL, SDDS_FLOAT) ||
      !SDDS_DefineSimpleColumn(SDDSout, "pr", NULL, SDDS_FLOAT) ||
      !SDDS_DefineSimpleColumn(SDDSout, "z0", "m", SDDS_FLOAT) ||
      !SDDS_DefineSimpleColumn(SDDSout, "r0", "m", SDDS_FLOAT) ||
      !SDDS_DefineSimpleColumn(SDDSout, "q", "C", SDDS_FLOAT) ||
      !SDDS_DefineSimpleColumn(SDDSout, "generation", NULL, SDDS_SHORT) ||
      !SDDS_DefineSimpleColumn(SDDSout, "particleID", NULL, SDDS_LONG) ||
#ifdef DEBUG
      !SDDS_DefineSimpleColumn(SDDSout, "lossCode", NULL, SDDS_SHORT) ||
      !SDDS_DefineSimpleColumn(SDDSout, "iz", NULL, SDDS_SHORT) ||
      !SDDS_DefineSimpleColumn(SDDSout, "ir", NULL, SDDS_SHORT) ||
      !SDDS_DefineSimpleColumn(SDDSout, "rLast", "m", SDDS_FLOAT) ||
      !SDDS_DefineSimpleColumn(SDDSout, "zLast", "m", SDDS_FLOAT) ||
#endif
      !SDDS_SaveLayout(SDDSout) || !SDDS_WriteLayout(SDDSout))
    SDDS_Bomb("problem setting up lost particle output file");
  if (filename != filename0)
    free(filename);
}

void setParticleVelocities(double *vz, double *vr, double *vphi, double *gamma,
                           double pz, double pr, double r_pphi, double r)
{
  if (r)
    {
      *gamma = sqrt(sqr(pz) + sqr(pr) + sqr(r_pphi / r) + 1);
      *vphi = c_mks * r_pphi / r / (*gamma);
    }
  else
    {
      *gamma = sqrt(sqr(pz) + sqr(pr) + 1);
      *vphi = 0;
    }
  *vz = c_mks * pz / (*gamma);
  *vr = c_mks * pr / (*gamma);
}

void interpolateFields(double *EzReturn, double *ErReturn, double *EphiReturn,
                       double *BzReturn, double *BrReturn, double *BphiReturn,
                       double z, double r, FIELDS *EM_problem)
{
  long iz, ir, nz, nr, iz1, ir1;
  long longit_interp, radial_interp;
  double dz, dr, fr, fl, fa, fb;
  double zmin;
  /*double zmax, rmax;*/
  double Ez, Er, Ephi, Bz, Br, Bphi;
  double EzOOE, ErOOE, BphiOOE;

  longit_interp = EM_problem->modes & FL_LONGIT_INTERP;
  radial_interp = EM_problem->modes & FL_RADIAL_INTERP;
  dz = EM_problem->dz;
  dr = EM_problem->dr;
  nz = EM_problem->nz;
  nr = EM_problem->nr;
  zmin = EM_problem->zmin;
  /*zmax = EM_problem->zmax;
    rmax = EM_problem->rmax;*/

  /* compute indices of closest grid point */
  iz = (z - zmin) / dz + 0.5;
  ir = r / dr + 0.5;

  Ez = Er = Bphi = Ephi = Br = Bz = 0;
  /* -- Ez and/or Br */
  if (longit_interp)
    {
      /* iz is to the left of the particle */
      iz = (fr = (z - zmin) / dz + 0.5);
      fr -= iz;
      fl = 1 - fr;
    }
  else
    {
      /* iz is closest to the particle */
      iz = (z - zmin) / dz + 1;
      fl = 1;
      fr = 0;
    }
  if (radial_interp)
    {
      /* ir is below the particle */
      ir = (fa = r / dr);
      fa -= ir;
      fb = 1 - fa;
    }
  else
    {
      /* ir is closest to the particle */
      ir = r / dr + 0.5;
      fa = 0;
      fb = 1;
    }
  if (iz >= 0 && iz < nz && ir >= 0 && ir < nr)
    {
      if ((iz1 = iz + 1) >= nz)
        iz1 = iz;
      if ((ir1 = ir + 1) >= nr)
        ir1 = ir;

      /* these lines take care of problems with Ez being zero on one side of a surface that runs
       * along a z=constant line that lies between two Ez mesh points.  This condition is indicated
       * by Er=0.  Could use quadratic interpolation for these values, but it probably isn't worth it.
       */
      if (EM_problem->metal_flags[iz][ir + 1] & FL_ER_ZERO && EM_problem->metal_flags[iz][ir] & FL_ER_ZERO)
        {
          if (iz > 0 && (EM_problem->metal_flags[iz - 1][ir] || EM_problem->metal_flags[iz - 1][ir + 1]))
            iz = iz1;
          else
            iz1 = iz;
        }
      else if (EM_problem->metal_flags[iz1][ir + 1] & FL_ER_ZERO && EM_problem->metal_flags[iz1][ir] & FL_ER_ZERO)
        {
          iz1 = iz;
        }
      if (EM_problem->modes & FL_TM_FIELDS)
        Ez = fb * (fl * (EM_problem->Ez[iz][ir] + EM_problem->EzImposed[iz][ir]) +
                   fr * (EM_problem->Ez[iz1][ir] + EM_problem->EzImposed[iz1][ir])) +
          fa * (fl * (EM_problem->Ez[iz][ir1] + EM_problem->EzImposed[iz][ir1]) +
                fr * (EM_problem->Ez[iz1][ir1] + EM_problem->EzImposed[iz1][ir1]));
#define ALL_FIELDS 1
#if defined(ALL_FIELDS)
      if (r != 0 && EM_problem->modes & FL_TE_FIELDS)
        Br = fb * (fl * EM_problem->Br2[iz][ir] + fr * EM_problem->Br2[iz1][ir]) + fa * (fl * EM_problem->Br2[iz][ir1] + fr * EM_problem->Br2[iz1][ir1]);
#endif
    }

#if defined(ALL_FIELDS)
  /* -- Er and/or Bz */
  if (longit_interp)
    {
      iz = (fr = (z - zmin) / dz);
      fr -= iz;
      fl = 1 - fr;
    }
  else
    {
      iz = (z - zmin) / dz + 0.5;
      fr = 0;
      fl = 1;
    }
  if (radial_interp)
    {
      ir = (fa = r / dr + 0.5);
      fa -= ir;
      fb = 1 - fa;
    }
  else
    {
      ir = r / dr + 1;
      fa = 0;
      fb = 1;
    }
  if (iz >= 0 && iz < nz && ir >= 0)
    {
      if (ir >= nr)
        ir = nr - 1;
      iz1 = iz + 1;
      ir1 = ir + 1;
      if (EM_problem->metal_flags[iz + 1][ir] & FL_EZ_ZERO && EM_problem->metal_flags[iz][ir] & FL_EZ_ZERO)
        ir = ir1;
      else if (EM_problem->metal_flags[iz + 1][ir + 1] & FL_EZ_ZERO && EM_problem->metal_flags[iz][ir + 1] & FL_EZ_ZERO)
        ir1 = ir;
      if (r != 0 && EM_problem->modes & FL_TM_FIELDS)
        Er = fl * (fb * (EM_problem->Er[iz][ir] + EM_problem->ErImposed[iz][ir]) + fa * (EM_problem->Er[iz][ir1] + EM_problem->ErImposed[iz][ir1])) + fr * (fb * (EM_problem->Er[iz1][ir] + EM_problem->ErImposed[iz1][ir]) + fa * (EM_problem->Er[iz1][ir1] + EM_problem->ErImposed[iz1][ir1]));
      if (EM_problem->modes & FL_TE_FIELDS)
        Bz = fl * (fb * EM_problem->Bz2[iz][ir] + fa * EM_problem->Bz2[iz][ir1]) + fr * (fb * EM_problem->Bz2[iz1][ir] + fa * EM_problem->Bz2[iz1][ir1]);
    }

  /* -- Ephi */
  if (longit_interp)
    {
      iz = (fr = (z - zmin) / dz);
      fr -= iz;
      fl = 1 - fr;
    }
  else
    {
      iz = (z - zmin) / dz + 0.5;
      fr = 0;
      fl = 1;
    }
  if (radial_interp)
    {
      ir = (fa = r / dr);
      fa -= ir;
      fb = 1 - fa;
    }
  else
    {
      ir = r / dr + 0.5;
      fa = 0;
      fb = 1;
    }
  if (r != 0 && EM_problem->modes & FL_TE_FIELDS && iz >= 0 && iz < nz && ir >= 0 && ir < nr)
    {
      if ((iz1 = iz + 1) >= nz)
        iz1 = iz;
      if ((ir1 = ir + 1) >= nr)
        ir1 = ir;
      Ephi = fl * fb * EM_problem->Ephi[iz][ir] + fr * fb * EM_problem->Ephi[iz1][ir] + fl * fa * EM_problem->Ephi[iz][ir1] + fr * fa * EM_problem->Ephi[iz1][ir1];
    }

  if (EM_problem->modes & FL_TM_FIELDS)
    {
      /* -- Bphi at time t */
      if (longit_interp)
        {
          iz = (fr = (z - zmin) / dz + 0.5);
          fr -= iz;
          fl = 1 - fr;
        }
      else
        {
          iz = (z - zmin) / dz + 1;
          fr = 0;
          fl = 1;
        }
      if (radial_interp)
        {
          ir = (fa = r / dr + 0.5);
          fa -= ir;
          fb = 1 - fa;
        }
      else
        {
          ir = r / dr + 1;
          fa = 0;
          fb = 1;
        }
      if (r != 0 && iz >= 0 && iz < nz)
        {
          if ((iz1 = iz + 1) >= nz)
            iz1 = iz;
          if (ir >= nr)
            ir = nr - 1;
          if ((ir1 = ir + 1) >= nr)
            ir1 = ir;
          Bphi = fl * fb * EM_problem->Bphi2[iz][ir] + fr * fb * EM_problem->Bphi2[iz1][ir] + fa * fl * EM_problem->Bphi2[iz][ir1] + fr * fa * EM_problem->Bphi2[iz1][ir1];
        }
    }
#endif /* ALL_FIELDS */

  /* add spatially-constant DC components */
  Ez += EM_problem->constantEz;
  Er += EM_problem->constantEr;
  Ephi += EM_problem->constantEphi;
  Bz += EM_problem->constantBz;
  Br += EM_problem->constantBr;
  Bphi += EM_problem->constantBphi;

  /* add fields from on-axis field profiles */
  EzOOE = ErOOE = BphiOOE = 0;
  addOffAxisExpansionFields(&EzOOE, &ErOOE, &BphiOOE,
                            &EM_problem->EzOOEMin, &EM_problem->EzOOEMax, z, r,
                            EM_problem->time,
                            EM_problem->onAxisField, EM_problem->onAxisFields, 0);
  Ez += EzOOE;
  Er += ErOOE;
  Bphi += BphiOOE;

  if (EM_problem->turnOffEz)
    Ez = 0;
  if (EM_problem->turnOffEr)
    Er = 0;
  if (EM_problem->turnOffEphi)
    Ephi = 0;
  if (EM_problem->turnOffBz)
    Bz = 0;
  if (EM_problem->turnOffBr)
    Br = 0;
  if (EM_problem->turnOffBphi)
    Bphi = 0;

  if (EM_problem->BzImposed)
    {
      /* add in solenoidal DC magnetic fields */
      iz = (fr = (z - zmin) / dz);
      ir = (fa = r / dr);
      if (ir >= 0 && iz >= 0 && ir < nr && iz < nz)
        {
          fa -= ir;
          fb = 1 - fa;
          fr -= iz;
          fl = 1 - fr;
          Bz += fl * (fb * EM_problem->BzImposed[iz][ir] + fa * EM_problem->BzImposed[iz][ir + 1]) +
            fr * (fb * EM_problem->BzImposed[iz + 1][ir] + fa * EM_problem->BzImposed[iz + 1][ir + 1]);
          Br += fl * (fb * EM_problem->BrImposed[iz][ir] + fa * EM_problem->BrImposed[iz][ir + 1]) +
            fr * (fb * EM_problem->BrImposed[iz + 1][ir] + fa * EM_problem->BrImposed[iz + 1][ir + 1]);
        }
    }

  if (EzReturn)
    *EzReturn = Ez;
  if (ErReturn)
    *ErReturn = Er;
  if (EphiReturn)
    *EphiReturn = Ephi;

  if (BzReturn)
    *BzReturn = Bz;
  if (BrReturn)
    *BrReturn = Br;
  if (BphiReturn)
    *BphiReturn = Bphi;
}
