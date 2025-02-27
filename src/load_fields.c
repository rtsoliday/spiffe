/*************************************************************************\
 * Copyright (c) 2002 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2002 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: load_fields.c
 * purpose: setup and execution of field reloading
 *
 * Michael Borland, 1992-1995
 */
#include "spiffe.h"

long checkLoadFile(SDDS_TABLE *inTable);
long checkLoadColumn(SDDS_TABLE *inTable, char *name, char *units);

void perform_field_loading(
                           FIELDS *EM_problem1,
                           NAMELIST_TEXT *nl_text)
{
#include "load_fields.h"
  FIELDS EM_problem0, *EM_problem;
  double Ez_max, z_max = 0, dz, dz1, dr, dr1, time;
  long nz = 0, nr = 0, iz, ir, nz1, nr1;
  long iz_start, iz_end, iz_offset, ir_end, readCode;
  SDDS_TABLE inTable;
  char TMfields, TEfields;

  EM_problem = &EM_problem0;
  EM_problem->metal_flags = NULL;
  EM_problem->Ez = EM_problem->Jz = EM_problem->Er = EM_problem->Jr =
    EM_problem->Bphi = EM_problem->Bphi2 = EM_problem->Ephi = EM_problem->Jphi =
    EM_problem->Br = EM_problem->Br2 = EM_problem->Bz = EM_problem->Bz2 = NULL;
  EM_problem->constantEz = EM_problem->constantBz = 0;
  EM_problem->zt = 0;

  filename = NULL;
  Ez_peak = time_threshold = 0;
  factor = 1;

  process_namelist(&load_fields, nl_text);
  print_namelist(stdout, &load_fields);

  if (filename == NULL)
    bomb("no filename given for field samples", NULL);
  if (Ez_peak == 0 && factor == 0)
    bomb("Ez_peak and factor cannot both be zero", NULL);

  if (!SDDS_InitializeInput(&inTable, filename))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
  if (!checkLoadFile(&inTable))
    bomb("load file is invalid--check presence, type, and units of data", NULL);

  while ((readCode = SDDS_ReadTable(&inTable)) > 0)
    {
      if (!SDDS_GetParameter(&inTable, "time", &time))
        {
          SDDS_SetError("unable to get value of parameter \"time\" for field loading");
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
        }
      if (time < time_threshold)
        continue;
      printf("Using fields from t = %es\n", time);
      if (!SDDS_GetParameter(&inTable, "version", &EM_problem->version))
        {
          SDDS_SetError("unable to get value of parameter \"version\" for field loading");
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
        }
      printf("File %s is version %ld\n", filename, EM_problem->version);
      if (!SDDS_GetParameters(&inTable,
                              "nr", &EM_problem->nr, "nz", &EM_problem->nz,
                              "dr", &EM_problem->dr, "dz", &EM_problem->dz,
                              "zmin", &EM_problem->zmin, "zmax", &EM_problem->zmax,
                              "rmax", &EM_problem->rmax, "zt", &EM_problem->zt,
                              "TMfields", &TMfields, "TEfields", &TEfields, NULL))
        {
          SDDS_SetError("unable to read parameters for field loading");
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
        }
      if (!(EM_problem->metal_flags = SDDS_GetMatrixFromColumn(&inTable, "metalFlags", EM_problem->nz + 1,
                                                               EM_problem->nr + 1, 0)))
        {
          SDDS_SetError("unable to get metalFlags column as matrix");
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
        }
      nz1 = (nz = EM_problem->nz) + 1;
      nr1 = (nr = EM_problem->nr) + 1;
      EM_problem->Psi = (double **)zarray_2d(sizeof(**EM_problem->Psi), nz1, nr1);
      EM_problem->Q = (double **)zarray_2d(sizeof(**EM_problem->Q), nz1, nr1);
      EM_problem->area_c = tmalloc(sizeof(*EM_problem->area_c) * nr1);
      EM_problem->area_o = tmalloc(sizeof(*EM_problem->area_o) * nr1);
      EM_problem->modes = 0;
      if (TMfields == 'y')
        {
          EM_problem->modes |= FL_TM_FIELDS;
          if (!(EM_problem->Ez = SDDS_GetDoubleMatrixFromColumn(&inTable, "Ez", nz1, nr1, 0)) ||
              !(EM_problem->Er = SDDS_GetDoubleMatrixFromColumn(&inTable, "Er", nz1, nr1, 0)) ||
              !(EM_problem->Jz = SDDS_GetDoubleMatrixFromColumn(&inTable, "Jz", nz1, nr1, 0)) ||
              !(EM_problem->Jr = SDDS_GetDoubleMatrixFromColumn(&inTable, "Jr", nz1, nr1, 0)) ||
              !(EM_problem->Bphi = SDDS_GetDoubleMatrixFromColumn(&inTable, "Bphi", nz1, nr1, 0)) ||
              !(EM_problem->Bphi2 = SDDS_GetDoubleMatrixFromColumn(&inTable, "Bphi2", nz1, nr1, 0)))
            {
              SDDS_SetError("unable to get TM field columns as matrices");
              SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
            }
        }
      if (TEfields == 'y')
        {
          EM_problem->modes |= FL_TE_FIELDS;
          if (!(EM_problem->Ephi = SDDS_GetDoubleMatrixFromColumn(&inTable, "Ephi", nz1, nr1, 0)) ||
              !(EM_problem->Jphi = SDDS_GetDoubleMatrixFromColumn(&inTable, "Jphi", nz1, nr1, 0)) ||
              !(EM_problem->Br = SDDS_GetDoubleMatrixFromColumn(&inTable, "Br", nz1, nr1, 0)) ||
              !(EM_problem->Br2 = SDDS_GetDoubleMatrixFromColumn(&inTable, "Br2", nz1, nr1, 0)) ||
              !(EM_problem->Bz = SDDS_GetDoubleMatrixFromColumn(&inTable, "Bz", nz1, nr1, 0)) ||
              !(EM_problem->Bz2 = SDDS_GetDoubleMatrixFromColumn(&inTable, "Bz2", nz1, nr1, 0)))
            {
              SDDS_SetError("unable to get TM field columns as matrices");
              SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
            }
        }
      break;
    }
  if (readCode < 1)
    bomb("unable to get data for field loading", NULL);
  if (!SDDS_Terminate(&inTable))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

  if (Ez_peak)
    {
      if (!(EM_problem->modes & FL_TM_FIELDS))
        bomb("can't scale Ez_peak because TM mode fields are not enabled", NULL);

      /* scan on-axis Ez and scale to desired value */
      Ez_max = 0;
      for (iz = 1; iz <= nz; iz++)
        {
          if (fabs(Ez_max) < fabs(EM_problem->Ez[iz][0]))
            {
              Ez_max = EM_problem->Ez[iz][0];
              z_max = (iz - 0.5) * EM_problem->dz + EM_problem->zmin;
            }
        }
      if (fabs(Ez_max) < 1e-10)
        bomb("extremal value of on-axis Ez is < 1e-10--can't scale up", NULL);
      printf("Extremal value of on-axis Ez is %e V at z = %e m.\n", Ez_max, z_max);
      printf("This will be scaled to %e V.\n", Ez_peak);
      factor = Ez_peak / Ez_max;
    }

  if (factor != 1)
    {
      /* scale all fields by specified factor */
      printf("Scaling all fields by factor %e .\n", factor);
      if (EM_problem->modes & FL_TM_FIELDS)
        for (iz = 0; iz <= nz; iz++)
          {
            for (ir = 0; ir <= nr; ir++)
              {
                EM_problem->Ez[iz][ir] *= factor;
                EM_problem->Jz[iz][ir] *= factor;
                EM_problem->Er[iz][ir] *= factor;
                EM_problem->Jr[iz][ir] *= factor;
                EM_problem->Bphi[iz][ir] *= factor;
                EM_problem->Bphi2[iz][ir] *= factor;
              }
          }
      if (EM_problem->modes & FL_TE_FIELDS)
        for (iz = 0; iz <= nz; iz++)
          {
            for (ir = 0; ir <= nr; ir++)
              {
                EM_problem->Ephi[iz][ir] *= factor;
                EM_problem->Jphi[iz][ir] *= factor;
                EM_problem->Br[iz][ir] *= factor;
                EM_problem->Br2[iz][ir] *= factor;
                EM_problem->Bz[iz][ir] *= factor;
                EM_problem->Bz2[iz][ir] *= factor;
              }
          }
    }

  /* add fields in to pre-existing data */
  iz_start = 0;
  iz_end = EM_problem1->nz;
  iz_offset = 0;
  ir_end = EM_problem1->nr;
  dz1 = (EM_problem1->zmax - EM_problem1->zmin) / (EM_problem1->nz - 1);
  dz = (EM_problem->zmax - EM_problem->zmin) / (EM_problem->nz - 1);
  dr1 = EM_problem1->rmax / (EM_problem1->nr - 1);
  dr = EM_problem->rmax / (EM_problem->nr - 1);
  if (EM_problem->nz != EM_problem1->nz || EM_problem->nr != EM_problem1->nr ||
      EM_problem->zmax != EM_problem1->zmax || EM_problem->zmin != EM_problem1->zmin ||
      EM_problem->rmax != EM_problem1->rmax)
    {
      printf("warning: possible geometry mismtach in overlay fields:\n");
      printf("existing data:  nz=%" PRId32 "  zmin=%.15e  zmax=%.15e  dz = %.15e\n", EM_problem1->nz, EM_problem1->zmin,
             EM_problem1->zmax, dz1);
      printf("                nr=%" PRId32 "  rmax=%.15e  dr = %.15e\n", EM_problem1->nr, EM_problem1->rmax, dr1);
      printf("overlay data :  nz=%" PRId32 "  zmin=%.15e  zmax=%.15e  dz = %.15e\n", EM_problem->nz, EM_problem->zmin,
             EM_problem->zmax, dz);
      printf("                nr=%" PRId32 "  rmax=%.15e  dr = %.15e\n", EM_problem->nr, EM_problem->rmax, dr);
      if (fabs(dz / dz1 - 1) > 1e-6)
        bomb("dz's not close enough", NULL);
      if (fabs(dr / dr1 - 1) > 1e-6)
        bomb("dr's not close enough", NULL);
      iz_offset = (EM_problem1->zmin - EM_problem->zmin) / dz + 0.5;
      if (fabs(iz_offset - (EM_problem1->zmin - EM_problem->zmin) / dz) > 1e-6)
        bomb("misalignment of boundary and grid lines from existing data and data in file", NULL);
      if (iz_offset < 0)
        iz_start = -iz_offset;
      if ((EM_problem1->nz + iz_offset) > EM_problem->nz)
        iz_end = EM_problem->nz - iz_offset;
      if (EM_problem1->nr < EM_problem->nr)
        ir_end = EM_problem->nr;
    }
  printf("overlay check:\nregion overlayed:  [%e, %e] x [0, %e]\ndata region:  [%e, %e] x [0, %e]\n",
         EM_problem1->zmin + iz_start * dz1, EM_problem1->zmin + (iz_end - 1) * dz1, (ir_end - 1) * dr1,
         EM_problem->zmin + (iz_start + iz_offset) * dz, EM_problem->zmin + (iz_end + iz_offset - 1) * dz, (ir_end - 1) * dr);
  if (EM_problem->modes & FL_TM_FIELDS)
    {
      if (!(EM_problem1->modes & FL_TM_FIELDS))
        bomb("attempt to overlay TM mode fields onto problem which does not contain those fields", NULL);
      for (iz = iz_start; iz <= iz_end; iz++)
        {
          for (ir = 0; ir <= ir_end; ir++)
            {
              EM_problem1->Ez[iz][ir] += EM_problem->Ez[iz + iz_offset][ir];
              EM_problem1->Jz[iz][ir] += EM_problem->Jz[iz + iz_offset][ir];
              EM_problem1->Er[iz][ir] += EM_problem->Er[iz + iz_offset][ir];
              EM_problem1->Jr[iz][ir] += EM_problem->Jr[iz + iz_offset][ir];
              EM_problem1->Bphi[iz][ir] += EM_problem->Bphi[iz + iz_offset][ir];
              EM_problem1->Bphi2[iz][ir] += EM_problem->Bphi2[iz + iz_offset][ir];
            }
        }
    }
  if (EM_problem->modes & FL_TE_FIELDS)
    {
      if (!(EM_problem->modes & FL_TE_FIELDS))
        bomb("attempt to overlay TE mode fields onto problem which does not contain those fields", NULL);
      for (iz = iz_start; iz <= iz_end; iz++)
        {
          for (ir = 0; ir <= ir_end; ir++)
            {
              EM_problem1->Ephi[iz][ir] += EM_problem->Ephi[iz + iz_offset][ir];
              EM_problem1->Jphi[iz][ir] += EM_problem->Jphi[iz + iz_offset][ir];
              EM_problem1->Br[iz][ir] += EM_problem->Br[iz + iz_offset][ir];
              EM_problem1->Br2[iz][ir] += EM_problem->Br2[iz + iz_offset][ir];
              EM_problem1->Bz[iz][ir] += EM_problem->Bz[iz + iz_offset][ir];
              EM_problem1->Bz2[iz][ir] += EM_problem->Bz2[iz + iz_offset][ir];
            }
        }
    }
  printf("fields added to pre-existing problem\n");
  fflush(stdout);
  /* free the arrays used for the loaded fields */
  free_zarray_2d((void **)EM_problem->metal_flags, nz + 1, nr + 1);
  free(EM_problem->area_c);
  free(EM_problem->area_o);
  if (EM_problem->modes & FL_TM_FIELDS)
    {
      free_zarray_2d((void **)EM_problem->Ez, nz + 1, nr + 1);
      free_zarray_2d((void **)EM_problem->Jz, nz + 1, nr + 1);
      free_zarray_2d((void **)EM_problem->Er, nz + 1, nr + 1);
      free_zarray_2d((void **)EM_problem->Jr, nz + 1, nr + 1);
      free_zarray_2d((void **)EM_problem->Bphi, nz + 1, nr + 1);
      free_zarray_2d((void **)EM_problem->Bphi2, nz + 1, nr + 1);
      free_zarray_2d((void **)EM_problem->Psi, nz + 1, nr + 1);
      free_zarray_2d((void **)EM_problem->Q, nz + 1, nr + 1);
    }
  if (EM_problem->modes & FL_TE_FIELDS)
    {
      free_zarray_2d((void **)EM_problem->Ephi, nz + 1, nr + 1);
      free_zarray_2d((void **)EM_problem->Jphi, nz + 1, nr + 1);
      free_zarray_2d((void **)EM_problem->Br, nz + 1, nr + 1);
      free_zarray_2d((void **)EM_problem->Br2, nz + 1, nr + 1);
      free_zarray_2d((void **)EM_problem->Bz, nz + 1, nr + 1);
      free_zarray_2d((void **)EM_problem->Bz2, nz + 1, nr + 1);
    }
}

long checkLoadFile(SDDS_TABLE *inTable)
{
  if (SDDS_CheckParameter(inTable, "time", "s", SDDS_DOUBLE, stderr) != SDDS_CHECK_OKAY ||
      SDDS_CheckParameter(inTable, "version", NULL, SDDS_LONG, stderr) != SDDS_CHECK_OKAY ||
      SDDS_CheckParameter(inTable, "nr", NULL, SDDS_LONG, stderr) != SDDS_CHECK_OKAY ||
      SDDS_CheckParameter(inTable, "nz", NULL, SDDS_LONG, stderr) != SDDS_CHECK_OKAY ||
      SDDS_CheckParameter(inTable, "dr", "m", SDDS_DOUBLE, stderr) != SDDS_CHECK_OKAY ||
      SDDS_CheckParameter(inTable, "dz", "m", SDDS_DOUBLE, stderr) != SDDS_CHECK_OKAY ||
      SDDS_CheckParameter(inTable, "zmin", "m", SDDS_DOUBLE, stderr) != SDDS_CHECK_OKAY ||
      SDDS_CheckParameter(inTable, "zmax", "m", SDDS_DOUBLE, stderr) != SDDS_CHECK_OKAY ||
      SDDS_CheckParameter(inTable, "rmax", "m", SDDS_DOUBLE, stderr) != SDDS_CHECK_OKAY ||
      SDDS_CheckParameter(inTable, "zt", "m", SDDS_DOUBLE, stderr) != SDDS_CHECK_OKAY ||
      SDDS_CheckParameter(inTable, "TMfields", NULL, SDDS_CHARACTER, stderr) != SDDS_CHECK_OKAY ||
      SDDS_CheckParameter(inTable, "TEfields", NULL, SDDS_CHARACTER, stderr) != SDDS_CHECK_OKAY ||
      SDDS_CheckColumn(inTable, "metalFlags", NULL, SDDS_LONG, stderr) != SDDS_CHECK_OKAY)
    return 0;
  if (!checkLoadColumn(inTable, "Ez", "V/m") || !checkLoadColumn(inTable, "Er", "V/m") ||
      !checkLoadColumn(inTable, "Jz", "A/m$a2$n") || !checkLoadColumn(inTable, "Jr", "A/m$a2$n") ||
      !checkLoadColumn(inTable, "Bphi", "T") || !checkLoadColumn(inTable, "Bphi2", "T"))
    return 0;
  if (!checkLoadColumn(inTable, "Ephi", "V/m") || !checkLoadColumn(inTable, "Jphi", "A/m$a2$n") ||
      !checkLoadColumn(inTable, "Br", "T") || !checkLoadColumn(inTable, "Br2", "T") ||
      !checkLoadColumn(inTable, "Bz", "T") || !checkLoadColumn(inTable, "Bz2", "T"))
    return 0;
  return 1;
}

long checkLoadColumn(SDDS_TABLE *inTable, char *name, char *units)
/* Check if the named column exists.  If it does, it must have the given units
 * and a floating point type.  It is okay if the column is nonexistent, however
 */
{
  long code;
  if ((code = SDDS_CheckColumn(inTable, name, units, SDDS_ANY_FLOATING_TYPE, NULL)) == SDDS_CHECK_OKAY)
    return 1;
  if (code == SDDS_CHECK_WRONGUNITS)
    return 0;
  fprintf(stdout, "Warning: column %s not found in field input file\n", name);
  return 1;
}
