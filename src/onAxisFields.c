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

void process_on_axis_field_definition(
                                      FIELDS *EM_problem,
                                      NAMELIST_TEXT *nl_text /* unparsed parameters in namelist text format */
                                      )
{
#include "onAxisFields.h"
  SDDS_DATASET SDDSin;
  long code, iz, ip;
  double Emin, Emax;

  if (!EM_problem)
    bomb("can't define on-axis fields before geometry is set up", NULL);

  /* read and check input parameters */
  filename = z_name = Ez_name = NULL;
  Ez_peak = frequency = phase = z_offset = 0;
  expansion_order = 3;

  process_namelist(&add_on_axis_fields, nl_text);
  print_namelist(stdout, &add_on_axis_fields);

  if (!filename || !strlen(filename))
    bomb("filename is blank", NULL);
  if (!z_name || !strlen(z_name))
    bomb("z_name is blank", NULL);
  if (!Ez_name || !strlen(Ez_name))
    bomb("Ez_name is blank", NULL);
  if (Ez_peak == 0)
    {
      fprintf(stdout, "*** Warning: Ez_peak=0, so this field is ignored.");
      return;
    }
  if (frequency < 0)
    bomb("frequency is negative", NULL);
  if (frequency == 0)
    fprintf(stdout, "*** Warning: frequency=0!");
  if (expansion_order < 0 || expansion_order > 3)
    bomb("expansion_order must by on [0, 3]", NULL);
  if (!(EM_problem->onAxisField = realloc(EM_problem->onAxisField,
                                          sizeof(*EM_problem->onAxisField) * (EM_problem->onAxisFields + 1))))
    bomb("memory allocation failure", NULL);

  if (!SDDS_InitializeInput(&SDDSin, filename))
    SDDS_Bomb("problem reading input file");

  if (SDDS_CheckColumn(&SDDSin, z_name, "m", SDDS_ANY_FLOATING_TYPE, stdout) != SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, Ez_name, NULL, SDDS_ANY_FLOATING_TYPE, stdout) != SDDS_CHECK_OK)
    SDDS_Bomb("invalid/missing columns in input file");

  /* read data from file */
  ip = EM_problem->onAxisFields;
  EM_problem->onAxisField[ip].frequency = frequency;
  EM_problem->onAxisField[ip].expansion_order = expansion_order;
  EM_problem->onAxisField[ip].phase = phase * PI / 180.0;
  while ((code = SDDS_ReadPage(&SDDSin)) > 0)
    {
      if (code > 1)
        {
          fprintf(stdout, "*** Warning: multiple pages in input file.  Only the first is used.\n");
          break;
        }
      if ((EM_problem->onAxisField[ip].points = SDDS_RowCount(&SDDSin)) <= 4)
        bomb("Too few points in file (need at least 4)", NULL);
      if (!(EM_problem->onAxisField[ip].z = SDDS_GetColumnInDoubles(&SDDSin, z_name)) ||
          !(EM_problem->onAxisField[ip].Ez = SDDS_GetColumnInDoubles(&SDDSin, Ez_name)))
        SDDS_Bomb("Problem getting data for on-axis field.");
      for (iz = 0; iz < EM_problem->onAxisField[ip].points; iz++)
        EM_problem->onAxisField[ip].z[iz] += z_offset;
      for (iz = 0; iz < EM_problem->onAxisField[ip].points - 1; iz++)
        {
          if (EM_problem->onAxisField[ip].z[iz + 1] <= EM_problem->onAxisField[ip].z[iz])
            bomb("Data points are not monotonically increasing in z.", NULL);
        }
    }
  if (code == 0)
    SDDS_Bomb("Problem with file");

  /* scale field to get desired peak value */
  find_min_max(&Emin, &Emax, EM_problem->onAxisField[ip].Ez, EM_problem->onAxisField[ip].points);
  Emin = fabs(Emin);
  Emax = fabs(Emax);
  if (Emax < Emin)
    Emax = Emin;
  if (Emax > 0)
    for (iz = 0; iz < EM_problem->onAxisField[ip].points; iz++)
      EM_problem->onAxisField[ip].Ez[iz] *= Ez_peak / Emax;

  /* take derivatives */
  if (!(EM_problem->onAxisField[ip].DzEz =
        malloc(sizeof(*(EM_problem->onAxisField[ip].DzEz)) * EM_problem->onAxisField[ip].points)) ||
      !(EM_problem->onAxisField[ip].Dz2Ez =
        malloc(sizeof(*(EM_problem->onAxisField[ip].Dz2Ez)) * EM_problem->onAxisField[ip].points)) ||
      !(EM_problem->onAxisField[ip].Dz3Ez =
        malloc(sizeof(*(EM_problem->onAxisField[ip].Dz3Ez)) * EM_problem->onAxisField[ip].points)))
    bomb("Memory allocation failure", NULL);

  take_derivative(EM_problem->onAxisField[ip].DzEz,
                  EM_problem->onAxisField[ip].Ez,
                  EM_problem->onAxisField[ip].z,
                  EM_problem->onAxisField[ip].points,
                  1, 1e-4);
  take_derivative(EM_problem->onAxisField[ip].Dz2Ez,
                  EM_problem->onAxisField[ip].DzEz,
                  EM_problem->onAxisField[ip].z,
                  EM_problem->onAxisField[ip].points,
                  1, 1e-4);
  take_derivative(EM_problem->onAxisField[ip].Dz3Ez,
                  EM_problem->onAxisField[ip].Dz2Ez,
                  EM_problem->onAxisField[ip].z,
                  EM_problem->onAxisField[ip].points,
                  1, 1e-4);

  EM_problem->onAxisFields++;

  if (fields_used && strlen(fields_used))
    {
      SDDS_DATASET SDDSout;
      long zIndex = 0, EzIndex = 0, ErIndex = 0, BphiIndex = 0;
      double Ez, Er, Bphi, dummy;
      char *filename;

      if (!SDDS_CopyString(&filename, fields_used) ||
          !(filename = compose_filename(filename, EM_problem->rootname)) ||
          !SDDS_InitializeOutput(&SDDSout, SDDS_BINARY, 0, NULL, NULL, filename) ||
          (zIndex = SDDS_DefineColumn(&SDDSout, "z", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)) < 0 ||
          (EzIndex = SDDS_DefineColumn(&SDDSout, "Ez", NULL, "V/m", NULL, NULL, SDDS_DOUBLE, 0)) < 0 ||
          (ErIndex = SDDS_DefineColumn(&SDDSout, "Er/r", NULL, "V/m$a2$n", NULL, NULL, SDDS_DOUBLE, 0)) < 0 ||
          (BphiIndex = SDDS_DefineColumn(&SDDSout, "Bphi/r", NULL, "T/m", NULL, NULL, SDDS_DOUBLE, 0)) < 0 ||
          !SDDS_WriteLayout(&SDDSout) ||
          !SDDS_StartPage(&SDDSout, EM_problem->nz))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
      free(filename);
      for (iz = 0; iz < EM_problem->nz; iz++)
        {
          Ez = Er = Bphi = dummy = 0;
          addOffAxisExpansionFields(&Ez, &Er, &Bphi, &dummy, &dummy,
                                    EM_problem->zmin + iz * EM_problem->dz, 1e-6, 0.0,
                                    EM_problem->onAxisField, EM_problem->onAxisFields, 1);
          if (!SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, iz,
                                 zIndex, EM_problem->zmin + iz * EM_problem->dz,
                                 EzIndex, Ez, ErIndex, Er / 1e-6, BphiIndex, Bphi / 1e-6, -1))
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
        }
      if (!SDDS_WritePage(&SDDSout) || !SDDS_Terminate(&SDDSout))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
    }
}

void take_derivative(double *dFdz, double *F, double *z, long n_pts, long check_dz, double dzFracLimit)
/* z values must be equi-spaced for this to be correct */
{
  long i;
  double dz, dz_min, dz_max;
  if (check_dz)
    {
      dz_min = dz_max = z[1] - z[0];
      for (i = 2; i < n_pts; i++)
        {
          if ((dz = z[i] - z[i - 1]) < dz_min)
            dz_min = dz;
          if (dz > dz_max)
            dz_max = dz;
        }
      if ((dz_max - dz_min) / ((dz_max + dz_min) / 2) > dzFracLimit)
        {
          fprintf(stdout, "Error: points are not equispaced with sufficient accuracy\n");
          fprintf(stdout, "dz_max = %e, dz_min = %e\n",
                  dz_max, dz_min);
          exit(1);
        }
    }

  for (i = 1; i < n_pts - 1; i++)
    dFdz[i] = (F[i + 1] - F[i - 1]) / (z[i + 1] - z[i - 1]);
  /* first and last point must be handled with two-point formula */
  dFdz[0] = (F[1] - F[0]) / (z[1] - z[0]);
  i = n_pts - 1;
  dFdz[i] = (F[i] - F[i - 1]) / (z[i] - z[i - 1]);
}

int compare_double(const void *d1, const void *d2)
{
  double diff;
  diff = (*((double *)d1)) - (*((double *)d2));
  if (diff < 0)
    return -1;
  if (diff > 0)
    return 1;
  return 0;
}

void addOffAxisExpansionFields(double *Ez, double *Er, double *Bphi,
                               double *EzMin, double *EzMax,
                               double z, double r, double t,
                               ON_AXIS_FIELD *onAxisField, long onAxisFields,
                               long profile)
{
  long ip, iz;
  double Ez0, DzEz, Dz2Ez, BF, k, dz, sin_wt, cos_wt, omega;

  for (ip = 0; ip < onAxisFields; ip++)
    {
      iz = binaryArraySearch((void *)onAxisField[ip].z, sizeof(*(onAxisField[ip].z)),
                             onAxisField[ip].points, (void *)&z, compare_double, 1);
      if (iz < 0 || iz >= (onAxisField[ip].points - 1))
        continue;
      if (!(onAxisField[ip].z[iz] <= z && onAxisField[ip].z[iz + 1] > z))
        {
          fprintf(stderr, "Problem with position search in addOffAxisExpansionFields\n");
          fprintf(stderr, "z = %e\n", z);
          fprintf(stderr, "iz = %ld, z[%ld] = %e, z[%ld] = %e\n",
                  iz, iz, onAxisField[ip].z[iz], iz + 1, onAxisField[ip].z[iz + 1]);
          continue;
        }
      omega = onAxisField[ip].frequency * PIx2;
      if (profile)
        {
          sin_wt = cos_wt = 1;
        }
      else
        {
          sin_wt = sin(omega * t + onAxisField[ip].phase);
          cos_wt = cos(omega * t + onAxisField[ip].phase);
        }
      BF = omega / sqr(c_mks) * cos_wt;
      k = omega / c_mks;
      dz = z - onAxisField[ip].z[iz];
      Ez0 = onAxisField[ip].Ez[iz] + onAxisField[ip].DzEz[iz] * dz;
      *Ez += Ez0 * sin_wt;
      *Er = 0;
      if (onAxisField[ip].expansion_order >= 1)
        {
          DzEz = onAxisField[ip].DzEz[iz] + onAxisField[ip].Dz2Ez[iz] * dz;
          *Er += -r / 2 * DzEz * sin_wt;
          *Bphi += r / 2 * Ez0 * BF;
          if (onAxisField[ip].expansion_order >= 2)
            {
              Dz2Ez = onAxisField[ip].Dz2Ez[iz] + onAxisField[ip].Dz3Ez[iz] * dz;
              *Ez += -ipow(r, 2) / 4. * (Dz2Ez + sqr(k) * Ez0) * sin_wt;
              if (onAxisField[ip].expansion_order >= 3)
                {
                  *Er += ipow(r, 3) / 16. * (onAxisField[ip].Dz3Ez[iz] + sqr(k) * DzEz) * sin_wt;
                  *Bphi += ipow(r, 3) / 16. * (Dz2Ez + sqr(k) * Ez0) * BF;
                }
            }
        }
      if (*EzMin > *Ez)
        *EzMin = *Ez;
      if (*EzMax < *Ez)
        *EzMax = *Ez;
    }
}
