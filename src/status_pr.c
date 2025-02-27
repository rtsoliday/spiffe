/*************************************************************************\
 * Copyright (c) 2002 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2002 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: status_pr.c
 * purpose: do status printouts
 *
 * Michael Borland, 1992 
 */
#include "spiffe.h"

static long oldNumParticles = 0;

void do_status_printouts(FIELDS *EM_problem, BEAM *beam, CATHODE *cathode,
                         long check_divergence, long space_charge,
                         SDDS_DATASET *SDDSout, long *outputRow)
{
  long iz, ir, ip;
  double min, max, value, sum, *fen;
  double Qerror, Qsurface, dQ, dr2, drdz;
  double Ar_inner, Ar_outer, Az;
  double KE, FE, TE;
  double zMax, rMax, zMin, rMin, pzMax, pzMin, prMin, prMax;
  double EzAbsMax, ErAbsMax;
  double Qtotal;
  double zAve, zRMS, zSD;
  double rAve, rRMS, rSD;
  double pzAve, pzRMS, pzSD;
  double prAve, prRMS, prSD;
  long secondaries;

  printf("status report for step %ld at time = %23.16e\n", EM_problem->time_step,
         EM_problem->time);
  KE = FE = TE = 0;
  zMax = rMax = zMin = rMin = pzMax = pzMin = prMin = prMax = 0;
  zAve = zRMS = zSD = 0;
  rAve = rRMS = rSD = 0;
  pzAve = pzRMS = pzSD = 0;
  prAve = prRMS = prSD = 0;

  EzAbsMax = ErAbsMax = Qtotal = 0;
  secondaries = 0;

  printf("\t%ld macroparticles in simulation, %ld of which are active\n", beam->np, beam->npActive);

  if (beam->np)
    {
      for (ip = 0; ip < beam->np; ip++)
        if (beam->status[ip] & PART_ACTIVE)
          {
            Qtotal += beam->Q[ip];
            if (beam->generation[ip] > 0)
              secondaries++;
          }

      printf("\tTotal active charge : %23.16e C\n", Qtotal);

      computeMoments(&zAve, &zRMS, &zSD, NULL, beam->z, beam->np);
      find_min_max(&min, &max, beam->z, beam->np);
      zMax = max;
      zMin = min;
      printf("\tExtremal    z: %23.16e, %23.16e m\n", min, max);

      computeMoments(&rAve, &rRMS, &rSD, NULL, beam->r, beam->np);
      find_min_max(&min, &max, beam->r, beam->np);
      rMin = min;
      rMax = max;
      printf("\tExtremal    r: %23.16e, %23.16e m\n", min, max);

      computeMoments(&pzAve, &pzRMS, &pzSD, NULL, beam->pz, beam->np);
      find_min_max(&min, &max, beam->pz, beam->np);
      pzMin = min;
      pzMax = max;
      printf("\tExtremal   Pz: %23.16e, %23.16e\n", min, max);

      computeMoments(&prAve, &prRMS, &prSD, NULL, beam->pr, beam->np);
      find_min_max(&min, &max, beam->pr, beam->np);
      prMin = min;
      prMax = max;
      printf("\tExtremal   Pr: %23.16e, %23.16e\n", min, max);

      find_min_max(&min, &max, beam->vz, beam->np);
      printf("\tExtremal   Vz: %23.16e, %23.16e m/s\n", min, max);

      find_min_max(&min, &max, beam->vr, beam->np);
      printf("\tExtremal   Vr: %23.16e, %23.16e m/s\n", min, max);

      KE = kinetic_energy(beam);
      printf("\tKinetic energy: %23.16e\n", KE);
    }

  if (EM_problem->modes & FL_TM_FIELDS)
    {
      if (!EM_problem->Ez)
        bomb("Ez is unexpectedly NULL in advance_field_solutions", NULL);
      if (!EM_problem->Jz)
        bomb("Jz is unexpectedly NULL in advance_field_solutions", NULL);
      if (!EM_problem->Er)
        bomb("Er is unexpectedly NULL in advance_field_solutions", NULL);
      if (!EM_problem->Jr)
        bomb("Jr is unexpectedly NULL in advance_field_solutions", NULL);
      if (!EM_problem->Bphi)
        bomb("Bphi is unexpectedly NULL in advance_field_solutions", NULL);
      if (!EM_problem->Bphi2)
        bomb("Bphi2 is unexpectedly NULL in advance_field_solutions", NULL);
      if (check_divergence)
        {
          Qerror = Qsurface = 0;
          dr2 = sqr(EM_problem->dr);
          drdz = EM_problem->dr * EM_problem->dz;
          for (ir = 0; ir < EM_problem->nr; ir++)
            {
              if (ir == 0)
                {
                  Ar_inner = Ar_outer = PI * drdz;
                  Az = dr2 / 4 * PI;
                }
              else
                {
                  Ar_inner = PI * (2 * ir - 1) * drdz;
                  Ar_outer = PI * (2 * ir + 1) * drdz;
                  Az = PI * dr2 / 4 * (isqr(2 * ir + 1) - isqr(2 * ir - 1));
                }
              for (iz = 0; iz < EM_problem->nz; iz++)
                {
                  if (ir == 0)
                    dQ = (EM_problem->Er[iz][ir + 1] * Ar_outer +
                          (EM_problem->Ez[iz + 1][ir] - EM_problem->Ez[iz][ir]) * Az) *
                      epsilon_o;
                  else
                    dQ = (EM_problem->Er[iz][ir + 1] * Ar_outer - EM_problem->Er[iz][ir] * Ar_inner +
                          (EM_problem->Ez[iz + 1][ir] - EM_problem->Ez[iz][ir]) * Az) *
                      epsilon_o;
                  if (EM_problem->metal_flags[iz][ir] & FL_IS_METAL || (iz == 0 && EM_problem->left_bc == NEUMANN) || (iz == EM_problem->nz - 1 && EM_problem->right_bc == NEUMANN) || (ir == EM_problem->nr - 1 && EM_problem->upper_bc == NEUMANN))
                    {
                      if (iz == 0)
                        dQ += EM_problem->Ez[iz][ir] * Az * epsilon_o;
                      if (iz == EM_problem->nz - 1)
                        dQ -= EM_problem->Ez[iz + 1][ir] * Az * epsilon_o;
                      if (ir == EM_problem->nr - 1)
                        dQ += EM_problem->Ez[iz][ir + 1] * Ar_outer * epsilon_o;
                      Qsurface += dQ;
                    }
                  else
                    Qerror += dQ;
                }
            }
          printf("\tError charge  : %23.16e C\n", Qerror - (space_charge && beam->np ? Qtotal : 0));
          printf("\tSurface charge: %23.16e C\n", Qsurface);
        }

      if (EM_problem->modes & FL_TE_FIELDS)
        printf("\t**** TM fields:\n");

      max = -(min = DBL_MAX);
      for (iz = 1; iz < EM_problem->nz; iz++)
        {
          for (ir = 0; ir < EM_problem->nr; ir++)
            {
              if ((value = EM_problem->Ez[iz][ir] + EM_problem->EzImposed[iz][ir]) > max)
                max = value;
              if (value < min)
                min = value;
            }
        }
      EzAbsMax = fabs(min) > fabs(max) ? fabs(min) : fabs(max);
      printf("\tExtremal   Ez: %23.16e, %23.16e V/m\n", min, max);

      max = -(min = DBL_MAX);
      sum = 0;
      for (iz = 1; iz < EM_problem->nz; iz++)
        {
          for (ir = 0; ir < EM_problem->nr; ir++)
            {
              if ((value = EM_problem->Jz[iz][ir]) > max)
                max = value;
              if (value < min)
                min = value;
              sum += value;
            }
        }
      printf("\tExtremal   Jz: %23.16e, %23.16e A/m^3\n", min, max);

      max = -(min = DBL_MAX);
      for (iz = 0; iz < EM_problem->nz; iz++)
        {
          for (ir = 0; ir < EM_problem->nr; ir++)
            {
              if ((value = EM_problem->Er[iz][ir] + EM_problem->ErImposed[iz][ir]) > max)
                max = value;
              if (value < min)
                min = value;
            }
        }
      ErAbsMax = fabs(min) > fabs(max) ? fabs(min) : fabs(max);
      printf("\tExtremal   Er: %23.16e, %23.16e V/m\n", min, max);

      max = -(min = DBL_MAX);
      sum = 0;
      for (iz = 0; iz < EM_problem->nz; iz++)
        {
          for (ir = 0; ir < EM_problem->nr; ir++)
            {
              if ((value = EM_problem->Jr[iz][ir]) > max)
                max = value;
              if (value < min)
                min = value;
              sum += value;
            }
        }
      printf("\tExtremal   Jr: %23.16e, %23.16e A/m^3\n", min, max);

      max = -(min = DBL_MAX);
      for (iz = 1; iz < EM_problem->nz; iz++)
        {
          for (ir = 1; ir < EM_problem->nr; ir++)
            {
              if ((value = EM_problem->Bphi[iz][ir]) > max)
                max = value;
              if (value < min)
                min = value;
            }
        }
      printf("\tExtremal Bphi: %23.16e, %23.16e T\n", min, max);
    }

  if (EM_problem->modes & FL_TE_FIELDS)
    {
      if (!EM_problem->Ephi)
        bomb("Ephi is unexpectedly NULL in advance_field_solutions", NULL);
      if (!EM_problem->Jphi)
        bomb("Jphi is unexpectedly NULL in advance_field_solutions", NULL);
      if (!EM_problem->Br)
        bomb("Br is unexpectedly NULL in advance_field_solutions", NULL);
      if (!EM_problem->Br2)
        bomb("Br2 is unexpectedly NULL in advance_field_solutions", NULL);
      if (!EM_problem->Bz)
        bomb("Bz is unexpectedly NULL in advance_field_solutions", NULL);
      if (!EM_problem->Bz2)
        bomb("Bz2 is unexpectedly NULL in advance_field_solutions", NULL);

      if (EM_problem->modes & FL_TM_FIELDS)
        printf("\t**** TE fields:\n");

      max = -(min = DBL_MAX);
      for (iz = 0; iz < EM_problem->nz; iz++)
        {
          for (ir = 0; ir < EM_problem->nr; ir++)
            {
              if ((value = EM_problem->Ephi[iz][ir]) > max)
                max = value;
              if (value < min)
                min = value;
            }
        }
      printf("\tExtremal Ephi: %23.16e, %23.16e V/m\n", min, max);

      max = -(min = DBL_MAX);
      for (iz = sum = 0; iz < EM_problem->nz; iz++)
        {
          for (ir = 0; ir < EM_problem->nr; ir++)
            {
              if ((value = EM_problem->Jphi[iz][ir]) > max)
                max = value;
              if (value < min)
                min = value;
              sum += EM_problem->Jphi[iz][ir] * EM_problem->dz *
                PI * sqr(EM_problem->dr) * (isqr(ir + 1) - isqr(ir));
            }
        }
      printf("\tExtremal Jphi: %23.16e, %23.16e A/m^3\n", min, max);
      printf("\tTotal    Iphi: %23.16e A\n", sum);

      max = -(min = DBL_MAX);
      for (iz = 0; iz < EM_problem->nz; iz++)
        {
          for (ir = 1; ir < EM_problem->nr; ir++)
            {
              if ((value = EM_problem->Bz[iz][ir]) > max)
                max = value;
              if (value < min)
                min = value;
            }
        }
      printf("\tExtremal   Bz: %23.16e, %23.16e T\n", min, max);

      max = -(min = DBL_MAX);
      for (iz = 1; iz < EM_problem->nz; iz++)
        {
          for (ir = 0; ir < EM_problem->nr; ir++)
            {
              if ((value = EM_problem->Br[iz][ir]) > max)
                max = value;
              if (value < min)
                min = value;
            }
        }
      printf("\tExtremal   Br: %23.16e, %23.16e T\n", min, max);
    }
  fen = field_energy(EM_problem);
  printf("\tField energy : %23.16e  (E: %23.16e, B: %23.16e)\n", FE = fen[0] + fen[1],
         fen[0], fen[1]);
  printf("\ttotal energy : %23.16e\n", TE = KE + FE);

  if (SDDSout && outputRow)
    {
      if (!SDDS_SetRowValues(SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, *outputRow,
                             "t", EM_problem->time,
                             "JCathode", (cathode ? cathode->current_density : 0),
                             "KineticEnergy", KE,
                             "FieldEnergy", FE,
                             "TotalEnergy", TE,
                             "Charge", Qtotal, "Particles", beam ? beam->np : 0,
                             "NewParticles", beam ? beam->np - oldNumParticles : 0,
                             "SecondaryParticles", secondaries,
                             "EzAbsMax", EzAbsMax, "ErAbsMax", ErAbsMax,
                             "zMax", zMax, "zMin", zMin,
                             "rMax", rMax, "rMin", rMin,
                             "pzMax", pzMax, "pzMin", pzMin,
                             "prMax", prMax, "prMin", prMin,
                             "zAve", zAve, "zRMS", zRMS, "zStDev", zSD,
                             "rAve", rAve, "rRMS", rRMS, "rStDev", rSD,
                             "pzAve", pzAve, "pzRMS", pzRMS, "pzStDev", pzSD,
                             "prAve", prAve, "prRMS", prRMS, "prStDev", prSD,
                             "EzOOEMin", EM_problem->EzOOEMin,
                             "EzOOEMax", EM_problem->EzOOEMax,
                             NULL) ||
          !SDDS_UpdatePage(SDDSout, FLUSH_TABLE))
        SDDS_Bomb("Problem writing status file");
      *outputRow += 1;
      oldNumParticles = beam->np;
    }

  fflush(stdout);
}

double kinetic_energy(BEAM *beam)
{
  long ip;
  double sum;
  double *pz, *pr, *r_pphi, *r;

  if (!beam->np)
    return (0.0);

  pz = beam->pz;
  pr = beam->pr;
  r = beam->r;
  r_pphi = beam->r_pphi;

  for (ip = sum = 0; ip < beam->np; ip++)
    sum += sqrt(1 + sqr(pz[ip]) + sqr(pr[ip]) + sqr(r[ip] ? r_pphi[ip] / r[ip] : 0));
  return (beam->stiffness * me_mks * sqr(c_mks) * (sum - beam->np));
}

double *field_energy(FIELDS *EM_problem)
{
  double **Ez, **Er, **Bphi, **Bz, **Br, **Ephi;
  double *area_o, *area_c;
  static double energy[2]; /* electric-field, magnetic-field energies */
  double sum, dz;
  /*double dr;*/
  long iz, ir, nr, nz;

  dz = EM_problem->dz;
  /*dr = EM_problem->dr;*/
  nr = EM_problem->nr;
  nz = EM_problem->nz;
  area_o = EM_problem->area_o;
  area_c = EM_problem->area_c;

  energy[0] = energy[1] = 0;
  if (EM_problem->modes & FL_TM_FIELDS)
    {
      Ez = EM_problem->Ez;
      Er = EM_problem->Er;
      Bphi = EM_problem->Bphi2;
      for (ir = 0; ir < nr; ir++)
        {
          sum = 0;
          for (iz = 1; iz < nz; iz++)
            {
              if (!(EM_problem->metal_flags[iz][ir] & FL_EZ_ZERO))
                sum += sqr(Ez[iz][ir]);
            }
          energy[0] += sum * area_c[ir] * dz * epsilon_o / 2;
        }
      for (ir = 1; ir < nr; ir++)
        {
          sum = 0;
          for (iz = 0; iz < nz; iz++)
            {
              if (!(EM_problem->metal_flags[iz][ir] & FL_ER_ZERO))
                sum += sqr(Er[iz][ir]);
            }
          energy[0] += sum * area_o[ir] * dz * epsilon_o / 2;
        }
      for (ir = 1; ir < nr; ir++)
        {
          sum = 0;
          for (iz = 1; iz < nz; iz++)
            {
              if (!(EM_problem->metal_flags[iz][ir] & FL_BPHI_ZERO))
                sum += sqr(Bphi[iz][ir]);
            }
          energy[1] += sum * area_o[ir] * dz / (2 * mu_o);
        }
    }

  if (EM_problem->modes & FL_TE_FIELDS)
    {
      Ephi = EM_problem->Ephi;
      Br = EM_problem->Br2;
      Bz = EM_problem->Bz2;
      for (ir = 0; ir < nr; ir++)
        {
          sum = 0;
          for (iz = 0; iz < nz; iz++)
            {
              if (!(EM_problem->metal_flags[iz][ir] & FL_EPHI_ZERO))
                sum += sqr(Ephi[iz][ir]);
            }
          energy[0] += sum * area_c[ir] * dz * epsilon_o / 2;
        }
      for (ir = 0; ir < nr; ir++)
        {
          sum = 0;
          for (iz = 1; iz < nz; iz++)
            {
              if (!(EM_problem->metal_flags[iz][ir] & FL_BR_ZERO))
                sum += sqr(Br[iz][ir]);
            }
          energy[1] += sum * area_c[ir] * dz / (2 * mu_o);
        }
      for (ir = 1; ir < nr; ir++)
        {
          sum = 0;
          for (iz = 0; iz < nz; iz++)
            {
              if (!(EM_problem->metal_flags[iz][ir] & FL_BZ_ZERO))
                sum += sqr(Bz[iz][ir]);
            }
          energy[1] += sum * area_o[ir] * dz / (2 * mu_o);
        }
    }
  return (energy);
}

void setUpStatusOutputFile(SDDS_DATASET *SDDSout, char *filename0, char *rootname)
{
  char *filename;
  filename = compose_filename(filename0, rootname);
  if (!SDDS_InitializeOutput(SDDSout, SDDS_BINARY, 0, NULL, NULL, filename) ||
      !SDDS_DefineSimpleColumn(SDDSout, "t", "s", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "JCathode", "A/m$a2$n", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "KineticEnergy", "J", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "FieldEnergy", "J", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "TotalEnergy", "J", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "Charge", "C", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "EzAbsMax", "V/m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "ErAbsMax", "V/m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "Particles", NULL, SDDS_LONG) ||
      !SDDS_DefineSimpleColumn(SDDSout, "NewParticles", NULL, SDDS_LONG) ||
      !SDDS_DefineSimpleColumn(SDDSout, "SecondaryParticles", NULL, SDDS_LONG) ||
      !SDDS_DefineSimpleColumn(SDDSout, "zMax", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "zMin", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "zAve", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "zRMS", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "zStDev", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "rMax", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "rMin", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "rAve", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "rRMS", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "rStDev", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "pzMax", NULL, SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "pzMin", NULL, SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "pzAve", NULL, SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "pzRMS", NULL, SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "pzStDev", NULL, SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "prMax", NULL, SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "prMin", NULL, SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "prAve", NULL, SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "prRMS", NULL, SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "prStDev", NULL, SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "EzOOEMin", "V/m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "EzOOEMax", "V/m", SDDS_DOUBLE) ||
      !SDDS_SaveLayout(SDDSout) || !SDDS_WriteLayout(SDDSout) ||
      !SDDS_StartPage(SDDSout, 1))
    SDDS_Bomb("problem setting up status output file");

  if (filename != filename0)
    free(filename);
  oldNumParticles = 0;
}
