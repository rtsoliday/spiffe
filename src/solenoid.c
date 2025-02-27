/*************************************************************************\
 * Copyright (c) 2002 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2002 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: solenoid.c
 * purpose: set up solenoidal fields.
 *
 * Michael Borland, 1995
 */
#include "spiffe.h"
#include "solenoid.h"

void process_solenoid_definition(
                                 FIELDS *EM_problem,
                                 NAMELIST_TEXT *nl_text /* unparsed parameters in namelist text format */
                                 )
{
  double **Bz, **Br, Bz1, Br1;
  long iz, ir, nz, nr, nrLimit;
  double z, r, dz, dr;
  static int32_t solenoidCount;

  if (!EM_problem)
    bomb("can't define a solenoid before geometry is set up", NULL);

  radius = evaluation_radius_limit = z_start = z_end = current = Bz_peak = 0;
  turns = 1;
  symmetry = 0;
  field_output = NULL;

  process_namelist(&define_solenoid, nl_text);
  print_namelist(stdout, &define_solenoid);
  if (radius <= 0)
    bomb("radius <= 0 is not valid", NULL);
  if (evaluation_radius_limit < 0)
    bomb("evaluation_radius_limit < 0  is not valid", NULL);
  if (z_start > z_end)
    bomb("z_start > z_end is not valid", NULL);
  if (turns <= 0)
    bomb("turns <= 0 is not valid", NULL);
  if (symmetry && symmetry != -1 && symmetry != 1)
    bomb("symmetry must be 0, 1, or -1", NULL);

  nz = EM_problem->nz;
  nr = EM_problem->nr;
  if (!(EM_problem->BzImposed))
    {
      EM_problem->BzImposed = (double **)zarray_2d(sizeof(**EM_problem->BzImposed), nz + 1, nr + 1);
      EM_problem->BrImposed = (double **)zarray_2d(sizeof(**EM_problem->BrImposed), nz + 1, nr + 1);
      solenoidCount = 0;
    }
  Bz = EM_problem->BzImposed;
  Br = EM_problem->BrImposed;

  if (!evaluation_radius_limit)
    nrLimit = nr;
  else if ((nrLimit = evaluation_radius_limit / EM_problem->dr + 0.5) > nr)
    nrLimit = nr;
  z = EM_problem->zmin;
  dz = EM_problem->dz;
  dr = EM_problem->dr;
  if (bucking)
    {
      iz = (z_buck - EM_problem->zmin) / EM_problem->dz;
      B_solenoid(&Br1, &Bz1, 0.0, z_buck, radius, z_start, z_end, turns, symmetry);
      if (Bz1 != 0)
        {
          current = -Bz[iz][0] / Bz1;
        }
      else
        {
          current = 0;
        }
      fprintf(stdout, "Bucking coil current set to %e\n", current);
    }
  for (iz = 0; iz <= nz; iz++)
    {
      r = 0;
      for (ir = 0; ir <= nrLimit; ir++)
        {
          B_solenoid(&Br1, &Bz1, r, z, radius, z_start, z_end, turns, symmetry);
          Br[iz][ir] += Br1 * current;
          Bz[iz][ir] += Bz1 * current;
          r += dr;
        }
      z += dz;
    }
  if (Bz_peak)
    {
      double Bz_max = -DBL_MAX;
      long iz_max = -1;
      for (iz = 0; iz <= nz; iz++)
        if (Bz[iz][0] > Bz_max)
          {
            iz_max = iz;
            Bz_max = Bz[iz][0];
          }
      if (iz_max == -1)
        bomb("problem scaling solenoid fields to Bz_peak---didn't find maximum!", NULL);
      else
        {
          if (Bz_max == 0)
            bomb("problem scaling solenoid fields to Bz_peak---maximum field is zero", NULL);
          fprintf(stdout, "Scaling fields by factor %e to achieve Bz_peak=%e T\n",
                  Bz_peak / Bz_max, Bz_peak);
          for (iz = 0; iz <= nz; iz++)
            {
              for (ir = 0; ir <= nrLimit; ir++)
                {
                  Bz[iz][ir] *= Bz_peak / Bz_max;
                  Br[iz][ir] *= Bz_peak / Bz_max;
                }
            }
        }
    }
  printf("solenoid computations finished\n");
  if (field_output)
    {
      dumpSolenoidFields(field_output, EM_problem->rootname, Br, Bz, EM_problem->zmin, dz, nz, dr, nr, ++solenoidCount);
      printf("accumulated solenoid field dump finished\n");
    }
}

void dumpSolenoidFields(char *filename0, char *rootname, double **Br, double **Bz, double zmin, double dz, long nz, double dr,
                        long nr, int32_t count)
{
  SDDS_TABLE outTable;
  long iz, ir, row;
  double r, z;
  long zIndex = 0, rIndex = 0, BrIndex = 0, BzIndex = 0;
  char *filename;

  if (!filename0)
    bomb("no filename passed to dumpSolenoidFields", NULL);
  if (!Br || !Bz)
    bomb("null field array passed to dumpSolenoidFields", NULL);
  if (nz <= 0 || nr <= 0)
    bomb("invalid dimensions passed to dumpSolenoidFields", NULL);
  filename = compose_filename(filename0, rootname);

  if (!SDDS_InitializeOutput(&outTable, SDDS_BINARY, 0, NULL, NULL, filename) ||
      (zIndex = SDDS_DefineColumn(&outTable, "z", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)) < 0 ||
      (rIndex = SDDS_DefineColumn(&outTable, "r", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)) < 0 ||
      (BzIndex = SDDS_DefineColumn(&outTable, "Bz", NULL, "T", NULL, NULL, SDDS_DOUBLE, 0)) < 0 ||
      (BrIndex = SDDS_DefineColumn(&outTable, "Br", NULL, "T", NULL, NULL, SDDS_DOUBLE, 0)) < 0 ||
      SDDS_DefineParameter1(&outTable, "SolenoidNumber", NULL, NULL, NULL, NULL, SDDS_LONG, &count) < 0 ||
      !SDDS_WriteLayout(&outTable) ||
      !SDDS_StartTable(&outTable, nz * nr))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
  if (filename != filename0)
    free(filename);

  z = zmin;
  for (iz = row = 0; iz < nz; iz++)
    {
      r = 0;
      for (ir = 0; ir < nr; ir++)
        {
          if (!SDDS_SetRowValues(&outTable, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE,
                                 row++,
                                 zIndex, z, rIndex, r,
                                 BzIndex, Bz[iz][ir], BrIndex, Br[iz][ir], -1))
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
          r += dr;
        }
      z += dz;
    }
  if (!SDDS_WriteTable(&outTable) || !SDDS_Terminate(&outTable))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
}
