/*************************************************************************\
 * Copyright (c) 2002 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2002 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: poisson_correction.c
 * purpose: do Poisson correction
 *
 * Michael Borland, 1992, 1994, 1995
 */
#include "spiffe.h"
#include "poisson_correction.h"

void poisson_guess(FIELDS *EM_problem, double charge, double centroid, long xguess_code);

#define NO_GUESS 0
#define LINE_CHARGE_GUESS 1
#define POINT_CHARGE_GUESS 2
#define ZERO_GUESS 3
#define N_GUESS_TYPES 4
static char *guess_name[N_GUESS_TYPES] = {
  "none", "line-charge", "point-charge", "zero"};

static long guess_code;
static long pc_start_step;

void setup_poisson_correction(
                              FIELDS *EM_problem,
                              NAMELIST_TEXT *nl_text /* unparsed parameters in namelist text format */
                              )
{

  /* set defaults */
  start_time = 0;
  step_interval = 0;
  accuracy = 1e-6;
  maximum_iterations = 1000;
  verbosity = 0;
  test_charge = z_test_charge = r_test_charge = 0;
  error_charge_threshold = 0;
  guess_type = "none";
  pc_start_step = -1;

  process_namelist(&poisson_correction, nl_text);
  print_namelist(stdout, &poisson_correction);

  if (start_time < 0)
    bomb("start_time < 0 for poisson_correction", NULL);
  if (step_interval < 0)
    bomb("step_interval < 0 for poisson_correction", NULL);
  if (accuracy <= 0)
    bomb("accuracy <= 0 for poisson_correction", NULL);
  if (maximum_iterations <= 0)
    bomb("maximum iterations <= 0 for poisson_correction", NULL);
  if (verbosity < 0)
    bomb("verbosity < 0 for poisson_correction", NULL);
  if (test_charge)
    {
      if (z_test_charge < EM_problem->zmin || z_test_charge > EM_problem->zmax ||
          r_test_charge < 0 || r_test_charge > EM_problem->rmax)
        bomb("test charge is outside of simulation region", NULL);
    }
  if ((guess_code = match_string(guess_type, guess_name, N_GUESS_TYPES, 0)) < 0)
    bomb("unknown guess type for poisson_correction", NULL);
  if (!(EM_problem->modes & FL_TM_FIELDS))
    bomb("TM fields must be enabled in order to do poisson correction", NULL);
}

void add_poisson_correction(FIELDS *EM_problem, BEAM *beam, long space_charge)
{
  long iz, ir;
  long ip, n_iter, longit_smear, radial_smear;
  double Az, Ar_inner, Ar_outer;
  double dr2, drdz, dr, dz, Qerror, Qreal, Qsurface, dQ;
  /*double dz2;*/
  double *volume;
  double centroid, fa, fb, fr, fl;
  BEAM testbeam;

  if (maximum_iterations <= 0)
    return;

  if (pc_start_step == -1)
    {
      if (EM_problem->dti <= 0)
        bomb("integration step size is not defined (add_poisson_correction)", NULL);
      pc_start_step = (start_time - EM_problem->start_time) / EM_problem->dti;
    }

  dr = EM_problem->dr;
  dz = EM_problem->dz;
  dr2 = sqr(dr);
  /*dz2 = sqr(dz);*/
  drdz = dr * dz;
  longit_smear = EM_problem->modes & FL_LONGIT_SMEAR;
  radial_smear = EM_problem->modes & FL_RADIAL_SMEAR;

  if (test_charge && EM_problem->time_step == 0)
    {
      printf("doing Poisson equation solution for test charge of %e C at z=%em, r=%em\n", test_charge,
             z_test_charge, r_test_charge);
      testbeam.max_np = testbeam.np = 1;
      testbeam.Q = tmalloc(sizeof(*testbeam.Q) * 1);
      testbeam.z = tmalloc(sizeof(*testbeam.z) * 1);
      testbeam.r = tmalloc(sizeof(*testbeam.r) * 1);
      testbeam.Q[0] = test_charge;
      testbeam.z[0] = z_test_charge;
      testbeam.r[0] = r_test_charge;
      beam = &testbeam;
    }

  if (step_interval <= 0 || EM_problem->time_step % step_interval != 0 || pc_start_step > EM_problem->time_step)
    {
      return;
    }

  if (!EM_problem->Psi)
    bomb("Psi array is unexpectedly null in add_poisson_correction", NULL);
  if (!EM_problem->Q)
    bomb("Q array is unexpectedly null in add_poisson_correction", NULL);
  if (!EM_problem->Ez)
    bomb("Ez array is unexpectedly null in add_poisson_correction", NULL);
  if (!EM_problem->Er)
    bomb("Er array is unexpectedly null in add_poisson_correction", NULL);

  if (verbosity)
    printf("Beginning Poisson correction\n");
  fflush(stdout);

  volume = tmalloc(sizeof(*volume) * (EM_problem->nr + 1));

  if (verbosity)
    printf("    computing apparent charge\n");
  fflush(stdout);
  Qreal = Qerror = Qsurface = centroid = 0;
  /* first compute the apparent charge density at each grid point using Ez and Er */
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
      volume[ir] = Az * dz;
      for (iz = 0; iz < EM_problem->nz; iz++)
        {
          EM_problem->Q[iz][ir] = 0;
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
            {
              EM_problem->Q[iz][ir] = dQ;
              centroid += dQ * (iz * dz + EM_problem->zmin);
            }
        }
    }

  if (space_charge)
    {
      if (verbosity)
        printf("    removing actual charge\n");
      fflush(stdout);
      /* remove the actual charge from each grid point */
      for (ip = 0; ip < beam->np; ip++)
        {
          if (!(beam->status[ip] & PART_ACTIVE))
            continue;
          if (longit_smear)
            iz = (beam->z[ip] - EM_problem->zmin) / dz;
          else
            iz = (beam->z[ip] - EM_problem->zmin) / dz + 0.5;
          if (radial_smear)
            ir = beam->r[ip] / dr;
          else
            ir = beam->r[ip] / dr + 0.5;
          if (iz < 0 || iz >= EM_problem->nz)
            continue;
          if (ir < 0 || ir >= EM_problem->nr)
            continue;
          /* fa and fb are the fractions contributed to the rows Above and Below */
          if (radial_smear && ir != EM_problem->nr - 1)
            fb = 1 - (fa = beam->r[ip] / dr - ir);
          else
            {
              fb = 1;
              fa = 0;
            }
          /* fl and fr are the fractions contributed to the rows to the Left and Right */
          if (longit_smear && iz != EM_problem->nz - 1)
            fl = 1 - (fr = (beam->z[ip] - EM_problem->zmin) / dz - iz);
          else
            {
              fl = 1;
              fr = 0;
            }
          if (!(EM_problem->metal_flags[iz][ir] & FL_IS_METAL || (iz == 0 && EM_problem->left_bc == NEUMANN) || (iz >= EM_problem->nz - 1 && EM_problem->right_bc == NEUMANN) || (ir >= EM_problem->nr - 1 && EM_problem->upper_bc == NEUMANN)))
            {
              EM_problem->Q[iz][ir] -= beam->Q[ip] * fb * fl;
              Qreal += beam->Q[ip] * fb * fl;
            }
          if (!(EM_problem->metal_flags[iz + 1][ir] & FL_IS_METAL || (iz == 0 && EM_problem->left_bc == NEUMANN) || (iz >= EM_problem->nz && EM_problem->right_bc == NEUMANN) || (ir >= EM_problem->nr - 1 && EM_problem->upper_bc == NEUMANN)))
            {
              EM_problem->Q[iz + 1][ir] -= beam->Q[ip] * fb * fr;
              Qreal += beam->Q[ip] * fb * fr;
            }
          if (!(EM_problem->metal_flags[iz][ir + 1] & FL_IS_METAL || (iz == 0 && EM_problem->left_bc == NEUMANN) || (iz >= EM_problem->nz - 1 && EM_problem->right_bc == NEUMANN) || (ir >= EM_problem->nr && EM_problem->upper_bc == NEUMANN)))
            {
              EM_problem->Q[iz][ir + 1] -= beam->Q[ip] * fa * fl;
              Qreal += beam->Q[ip] * fa * fl;
            }
          if (!(EM_problem->metal_flags[iz + 1][ir + 1] & FL_IS_METAL || (iz == 0 && EM_problem->left_bc == NEUMANN) || (iz >= EM_problem->nz && EM_problem->right_bc == NEUMANN) || (ir >= EM_problem->nr && EM_problem->upper_bc == NEUMANN)))
            {
              EM_problem->Q[iz + 1][ir + 1] -= beam->Q[ip] * fa * fr;
              Qreal += beam->Q[ip] * fa * fr;
            }
        }
    }
  else
    Qreal = 0;

  Qerror = 0;
  for (iz = 0; iz < EM_problem->nz; iz++)
    for (ir = 0; ir < EM_problem->nr; ir++)
      Qerror += EM_problem->Q[iz][ir];

  if (verbosity)
    {
      printf("    \"real\" charge: %23.16e C    error charge: %23.16e C    surface charge: %23.16e\n",
             Qreal, Qerror, Qsurface);
    }

  if ((error_charge_threshold > 0 && fabs(Qerror) < error_charge_threshold) ||
      (error_charge_threshold < 0 && (Qreal == 0 || fabs(Qerror / Qreal) < fabs(error_charge_threshold))))
    {
      if (verbosity)
        puts("    error charge is below threshold--Poisson correction not done.");
      free(volume);
      return;
    }

  /* go through and convert Q to -rho/epsilon_o */
  for (ir = 0; ir < EM_problem->nr; ir++)
    {
      for (iz = 0; iz < EM_problem->nz; iz++)
        {
          EM_problem->Q[iz][ir] = -EM_problem->Q[iz][ir] / (volume[ir] * epsilon_o);
        }
    }

  /* solve Poisson's equation for the error potential */
  poisson_guess(EM_problem, Qerror, centroid, guess_code);
  n_iter = solve_poisson_cyl(EM_problem->Psi, EM_problem->Q,
                             EM_problem->metal_flags, EM_problem->left_bc, EM_problem->right_bc, EM_problem->nz, EM_problem->nr,
                             dz, dr, accuracy, maximum_iterations);
  if (verbosity)
    {
      fflush(stdout);
      printf("    Poisson solution took %ld iterations\n", n_iter);
      fflush(stdout);
    }

  /* Go through and subtract off the fields from the error potential. */
  /* -- Ez */
  for (ir = 0; ir < EM_problem->nr; ir++)
    {
      for (iz = 1; iz < EM_problem->nz; iz++)
        {
          if (!(EM_problem->metal_flags[iz][ir] & FL_EZ_ZERO))
            EM_problem->Ez[iz][ir] -= -(EM_problem->Psi[iz][ir] - EM_problem->Psi[iz - 1][ir]) / dz;
        }
      if (EM_problem->left_bc == NEUMANN)
        EM_problem->Ez[0][ir] = EM_problem->Ez[1][ir];
      else
        EM_problem->Ez[0][ir] = -EM_problem->Ez[1][ir];
      if (EM_problem->right_bc == NEUMANN)
        EM_problem->Ez[EM_problem->nz][ir] = EM_problem->Ez[EM_problem->nz - 1][ir];
      else
        EM_problem->Ez[EM_problem->nz][ir] = -EM_problem->Ez[EM_problem->nz - 1][ir];
    }
  /* Neumann B.C. on upper boundary */
  for (iz = 0; iz <= EM_problem->nz; iz++)
    EM_problem->Ez[iz][EM_problem->nr - 1] = 0;

  /* -- Er */
  for (ir = 1; ir <= EM_problem->nr; ir++)
    {
      for (iz = 0; iz < EM_problem->nz; iz++)
        {
          if (!(EM_problem->metal_flags[iz][ir] & FL_ER_ZERO))
            EM_problem->Er[iz][ir] -= -(EM_problem->Psi[iz][ir] - EM_problem->Psi[iz][ir - 1]) / dr;
        }
    }
  for (iz = 0; iz < EM_problem->nz; iz++)
    EM_problem->Er[iz][0] = -EM_problem->Er[iz][1]; /* Dirichlet conditions on lower boundary */
  if (EM_problem->left_bc == NEUMANN)
    for (ir = 0; ir <= EM_problem->nr; ir++)
      EM_problem->Er[0][ir] = 0;
  if (EM_problem->right_bc == NEUMANN)
    for (ir = 0; ir <= EM_problem->nr; ir++)
      EM_problem->Er[EM_problem->nz - 1][ir] = 0;

  if (verbosity)
    printf("    re-computing apparent charge\n");
  fflush(stdout);

  /* recompute the Q array (is this necessary? can I just save it from above?) */
  Qreal = Qerror = Qsurface = 0;
  /* first compute the apparent charge density at each grid point using Ez and Er */
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
      volume[ir] = Az * dz;
      for (iz = 0; iz < EM_problem->nz; iz++)
        {
          EM_problem->Q[iz][ir] = 0;
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
                dQ -= EM_problem->Ez[iz][ir + 1] * Ar_outer * epsilon_o;
              Qsurface += dQ;
            }
          else
            {
              EM_problem->Q[iz][ir] = dQ;
            }
        }
    }

  if (space_charge)
    {
      if (verbosity)
        printf("    removing actual charge\n");
      fflush(stdout);
      /* remove the actual charge from each grid point */
      for (ip = 0; ip < beam->np; ip++)
        {
          if (!(beam->status[ip] & PART_ACTIVE))
            continue;
          if (longit_smear)
            iz = (beam->z[ip] - EM_problem->zmin) / dz;
          else
            iz = (beam->z[ip] - EM_problem->zmin) / dz + 0.5;
          if (radial_smear)
            ir = beam->r[ip] / dr;
          else
            ir = beam->r[ip] / dr + 0.5;

          if (iz < 0 || iz >= EM_problem->nz)
            continue;
          if (ir < 0 || ir >= EM_problem->nr)
            continue;
          /* fa and fb are the fractions contributed to the rows Above and Below */
          if (radial_smear && ir != EM_problem->nr - 1)
            fb = 1 - (fa = beam->r[ip] / dr - ir);
          else
            {
              fb = 1;
              fa = 0;
            }
          /* fl and fr are the fractions contributed to the rows to the Left and Right */
          if (longit_smear && iz != EM_problem->nz - 1)
            fl = 1 - (fr = (beam->z[ip] - EM_problem->zmin) / dz - iz);
          else
            {
              fl = 1;
              fr = 0;
            }
          if (!(EM_problem->metal_flags[iz][ir] & FL_IS_METAL || (iz == 0 && EM_problem->left_bc == NEUMANN) || (iz >= EM_problem->nz - 1 && EM_problem->right_bc == NEUMANN) || (ir >= EM_problem->nr - 1 && EM_problem->upper_bc == NEUMANN)))
            {
              EM_problem->Q[iz][ir] -= beam->Q[ip] * fb * fl;
              Qreal += beam->Q[ip] * fb * fl;
            }
          if (!(EM_problem->metal_flags[iz + 1][ir] & FL_IS_METAL || (iz == 0 && EM_problem->left_bc == NEUMANN) || (iz >= EM_problem->nz && EM_problem->right_bc == NEUMANN) || (ir >= EM_problem->nr - 1 && EM_problem->upper_bc == NEUMANN)))
            {
              EM_problem->Q[iz + 1][ir] -= beam->Q[ip] * fb * fr;
              Qreal += beam->Q[ip] * fb * fr;
            }
          if (!(EM_problem->metal_flags[iz][ir + 1] & FL_IS_METAL || (iz == 0 && EM_problem->left_bc == NEUMANN) || (iz >= EM_problem->nz - 1 && EM_problem->right_bc == NEUMANN) || (ir >= EM_problem->nr && EM_problem->upper_bc == NEUMANN)))
            {
              EM_problem->Q[iz][ir + 1] -= beam->Q[ip] * fa * fl;
              Qreal += beam->Q[ip] * fa * fl;
            }
          if (!(EM_problem->metal_flags[iz + 1][ir + 1] & FL_IS_METAL || (iz == 0 && EM_problem->left_bc == NEUMANN) || (iz >= EM_problem->nz && EM_problem->right_bc == NEUMANN) || (ir >= EM_problem->nr && EM_problem->upper_bc == NEUMANN)))
            {
              EM_problem->Q[iz + 1][ir + 1] -= beam->Q[ip] * fa * fr;
              Qreal += beam->Q[ip] * fa * fr;
            }
        }
    }
  else
    Qreal = 0;

  if (verbosity)
    {
      Qerror = 0;
      for (iz = 0; iz < EM_problem->nz; iz++)
        for (ir = 0; ir < EM_problem->nr; ir++)
          Qerror += EM_problem->Q[iz][ir];

      printf("    \"real\" charge: %23.16e C    error charge: %23.16e C    surface charge: %23.16e\n", Qreal, Qerror, Qsurface);
      printf("    Poisson correction completed\n");
      fflush(stdout);
    }

  free(volume);
}

void compute_imposed_Efield(FIELDS *EM_problem, long initialize)
{
  long iz, ir;
  long n_iter;
  double dr, dz, field_factor;
  /*long longit_smear, radial_smear;
    double dr2, dz2, drdz;*/
  static double last_field_factor = -1;

  if (maximum_iterations <= 0)
    return;

  if (pc_start_step == -1)
    {
      if (EM_problem->dti <= 0)
        bomb("integration step size is not defined (compute_imposed_Efield)", NULL);
      pc_start_step = (start_time - EM_problem->start_time) / EM_problem->dti;
    }

  dr = EM_problem->dr;
  dz = EM_problem->dz;
  /*dr2 = sqr(dr);
    dz2 = sqr(dz);
    drdz = dr*dz;
    longit_smear = EM_problem->modes&FL_LONGIT_SMEAR;
    radial_smear = EM_problem->modes&FL_RADIAL_SMEAR;*/

  if (!EM_problem->PsiImposed)
    bomb("PsiImposed array is unexpectedly null in compute_imposed_Efield", NULL);
  if (!EM_problem->EzImposed)
    bomb("EzImposed array is unexpectedly null in compute_imposed_Efield", NULL);
  if (!EM_problem->ErImposed)
    bomb("ErImposed array is unexpectedly null in compute_imposed_Efield", NULL);

  if (initialize)
    {
      /* solve Poisson's equation for the imposed potential */
      n_iter = solve_poisson_cyl(EM_problem->PsiImposed, NULL,
                                 EM_problem->metal_flags, EM_problem->left_bc, EM_problem->right_bc, EM_problem->nz, EM_problem->nr,
                                 dz, dr, accuracy, maximum_iterations);
      if (verbosity)
        {
          fflush(stdout);
          printf("    Poisson solution for imposed potential took %ld iterations\n", n_iter);
          fflush(stdout);
        }
    }

  field_factor = 1;
  if (EM_problem->imposedEFieldRampTime)
    {
      if (EM_problem->time < EM_problem->imposedEFieldRampTime)
        field_factor = EM_problem->time / EM_problem->imposedEFieldRampTime;
      else if (EM_problem->time < (EM_problem->imposedEFieldRampTime + EM_problem->imposedEFieldFlatTopTime))
        field_factor = 1;
      else if (EM_problem->time < (2 * EM_problem->imposedEFieldRampTime + EM_problem->imposedEFieldFlatTopTime))
        field_factor = ((2 * EM_problem->imposedEFieldRampTime + EM_problem->imposedEFieldFlatTopTime) - EM_problem->time) / EM_problem->imposedEFieldRampTime;
      else
        field_factor = 0;
      if (field_factor < 0)
        field_factor = 0;
    }

  if (field_factor == last_field_factor)
    return;
  last_field_factor = field_factor;

  /* Go through and compute the fields from the potential */
  /* -- Ez */
  for (ir = 0; ir < EM_problem->nr; ir++)
    {
      for (iz = 1; iz < EM_problem->nz; iz++)
        {
          if (!(EM_problem->metal_flags[iz][ir] & FL_EZ_ZERO))
            EM_problem->EzImposed[iz][ir] = -field_factor * (EM_problem->PsiImposed[iz][ir] - EM_problem->PsiImposed[iz - 1][ir]) / dz;
        }
      if (EM_problem->left_bc == NEUMANN)
        EM_problem->EzImposed[0][ir] = EM_problem->EzImposed[1][ir];
      else
        EM_problem->EzImposed[0][ir] = -EM_problem->EzImposed[1][ir];
      if (EM_problem->right_bc == NEUMANN)
        EM_problem->EzImposed[EM_problem->nz][ir] = EM_problem->EzImposed[EM_problem->nz - 1][ir];
      else
        EM_problem->EzImposed[EM_problem->nz][ir] = -EM_problem->EzImposed[EM_problem->nz - 1][ir];
    }
  /* Neumann B.C. on upper boundary */
  for (iz = 0; iz <= EM_problem->nz; iz++)
    EM_problem->EzImposed[iz][EM_problem->nr - 1] = 0;

  /* -- Er */
  for (ir = 1; ir <= EM_problem->nr; ir++)
    {
      for (iz = 0; iz < EM_problem->nz; iz++)
        {
          if (!(EM_problem->metal_flags[iz][ir] & FL_ER_ZERO))
            EM_problem->ErImposed[iz][ir] = -field_factor * (EM_problem->PsiImposed[iz][ir] - EM_problem->PsiImposed[iz][ir - 1]) / dr;
        }
    }
  for (iz = 0; iz < EM_problem->nz; iz++)
    EM_problem->ErImposed[iz][0] = -EM_problem->ErImposed[iz][1]; /* Dirichlet conditions on lower boundary */
  if (EM_problem->left_bc == NEUMANN)
    for (ir = 0; ir <= EM_problem->nr; ir++)
      EM_problem->ErImposed[0][ir] = 0;
  if (EM_problem->right_bc == NEUMANN)
    for (ir = 0; ir <= EM_problem->nr; ir++)
      EM_problem->ErImposed[EM_problem->nz - 1][ir] = 0;
}

void poisson_guess(
                   FIELDS *EM_problem,
                   double charge,   /* total charge in problem */
                   double centroid, /* z coordinate of centroid of charge */
                   long xguess_code)
{
  long iz, ir, nz, nr;
  double dist, guess, dz, dr, zmin;
  double **Psi;

  zmin = EM_problem->zmin;
  dz = EM_problem->dz;
  dr = EM_problem->dr;
  nz = EM_problem->nz;
  nr = EM_problem->nr;
  Psi = EM_problem->Psi;

  switch (xguess_code)
    {
    case LINE_CHARGE_GUESS:
      /* put in Phi = Q/(2*PI*eps)*ln(r) as guess */
      for (ir = 0; ir < nr; ir++)
        {
          guess = charge / (PIx2 * epsilon_o) * log((ir ? ir : 0.5) * dr);
          for (iz = 0; iz < nz; iz++)
            Psi[iz][ir] = guess;
        }
      break;
    case POINT_CHARGE_GUESS:
      /* put in Phi = Q/(4*PI*eps*r) as guess */
      for (ir = 0; ir < nr; ir++)
        {
          for (iz = 0; iz < nz; iz++)
            {
              dist = sqrt(sqr((iz * dz + zmin) - centroid) + sqr(dr * (ir + 0.5)));
              Psi[iz][ir] = charge / (4 * PI * epsilon_o * dist);
            }
        }
      break;
    case ZERO_GUESS:
      for (iz = 0; iz < nz; iz++)
        {
          for (ir = 0; ir < nr; ir++)
            {
              Psi[iz][ir] = 0;
            }
        }
      break;
    default:
      break;
    }
}
