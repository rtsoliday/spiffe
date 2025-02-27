/*************************************************************************\
 * Copyright (c) 2002 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2002 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: field_saving.c
 * purpose: setup and execution of field saving
 *
 * Michael Borland, 1994-1995
 */
#include "spiffe.h"

void process_field_output_definition(
                                     FIELD_OUTPUT *output_defn,
                                     FIELDS *EM_problem,
                                     NAMELIST_TEXT *nl_text)
{
#include "field_output.h"

  filename = NULL;
  z_interval = r_interval = 1;
  time_interval = start_time = 0;
  exclude_imposed_field = separate_imposed_field = 0;

  process_namelist(&define_field_output, nl_text);
  print_namelist(stdout, &define_field_output);

  if (filename == NULL)
    bomb("no filename given for field samples", NULL);
  output_defn->filename = filename;
  if (time_interval < 0)
    bomb("time_interval < 0 for field samples", NULL);
  output_defn->time_interval = time_interval;
  output_defn->start_time = start_time;
  if ((output_defn->z_interval = z_interval) < 1 || (output_defn->r_interval = r_interval) < 1)
    bomb("z_interval and r_interval must be > 1", NULL);

  /* these will be set on the first call to integrate(), when the integration time-step is known */
  output_defn->step_interval = output_defn->start_step = -1;
  output_defn->output_index = -1;
  output_defn->exclude_imposed_field = exclude_imposed_field;
  output_defn->separate_imposed_field = separate_imposed_field;
}

void output_fields(
                   FIELD_OUTPUT *output_defn,
                   FIELDS *EM_problem)
{
  long iz, ir, row;
  long iEz, iEr, iBphi, iEphi, iBr, iBz, rIndex, zIndex;
  long iJz, iJr, iJphi, iRho, iEzImposed, iErImposed;
  float rValue, zValue, EzValue, ErValue, BphiValue, EphiValue, BrValue, BzValue;
  float JzValue, JrValue, JphiValue, EzImposedValue, ErImposedValue;
  char *filename;

  if (output_defn->filename == NULL)
    return;

  if (output_defn->output_index == -1)
    {
      output_defn->step_interval = output_defn->time_interval / EM_problem->dti + 0.5;
      output_defn->start_step = (output_defn->start_time - EM_problem->start_time) / EM_problem->dti + 0.5;
      output_defn->output_index = 0;
      filename = compose_filename(output_defn->filename, EM_problem->rootname);
      if (!SDDS_InitializeOutput(&output_defn->SDDSout, SDDS_BINARY, 0, NULL, "spiffe field output",
                                 filename) ||
          SDDS_DefineParameter(&output_defn->SDDSout, "Variable1Name", NULL, NULL, NULL, NULL, SDDS_STRING, "z") < 0 ||
          SDDS_DefineParameter(&output_defn->SDDSout, "zMinimum", NULL, "m", NULL, NULL, SDDS_DOUBLE, NULL) < 0 ||
          SDDS_DefineParameter(&output_defn->SDDSout, "zInterval", NULL, "m", NULL, NULL, SDDS_DOUBLE, NULL) < 0 ||
          SDDS_DefineParameter(&output_defn->SDDSout, "zDimension", NULL, "m", NULL, NULL, SDDS_LONG, NULL) < 0 ||
          SDDS_DefineParameter(&output_defn->SDDSout, "rMinimum", NULL, "m", NULL, NULL, SDDS_DOUBLE, NULL) < 0 ||
          SDDS_DefineParameter(&output_defn->SDDSout, "rInterval", NULL, "m", NULL, NULL, SDDS_DOUBLE, NULL) < 0 ||
          SDDS_DefineParameter(&output_defn->SDDSout, "rDimension", NULL, "m", NULL, NULL, SDDS_LONG, NULL) < 0 ||
          SDDS_DefineParameter(&output_defn->SDDSout, "Variable2Name", NULL, NULL, NULL, NULL, SDDS_STRING, "r") < 0 ||
          SDDS_DefineParameter(&output_defn->SDDSout, "Time", NULL, "s", "Time of output", NULL, SDDS_DOUBLE, NULL) < 0 ||
          SDDS_DefineParameter(&output_defn->SDDSout, "Timestep", NULL, NULL, "Time-step of output", NULL, SDDS_LONG, NULL) < 0)
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
      if (filename != EM_problem->rootname)
        free(filename);
      if ((zIndex = SDDS_DefineColumn(&output_defn->SDDSout, "z", NULL, "m", NULL, NULL, SDDS_FLOAT, 0)) < 0 ||
          (rIndex = SDDS_DefineColumn(&output_defn->SDDSout, "r", NULL, "m", NULL, NULL, SDDS_FLOAT, 0)) < 0 ||
          (iRho = SDDS_DefineColumn(&output_defn->SDDSout, "rho", "$gr$r", "C/m$a3$n", NULL, NULL, SDDS_FLOAT, 0)) < 0)
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
      if (EM_problem->modes & FL_TM_FIELDS)
        {
          if ((iEz = SDDS_DefineColumn(&output_defn->SDDSout, "Ez", "E$bz$n", "V/m", NULL, NULL, SDDS_FLOAT, 0)) < 0 ||
              (iEr = SDDS_DefineColumn(&output_defn->SDDSout, "Er", "E$br$n", "V/m", NULL, NULL, SDDS_FLOAT, 0)) < 0 ||
              (iBphi = SDDS_DefineColumn(&output_defn->SDDSout, "Bphi", "B$bF$n", "T", NULL, NULL, SDDS_FLOAT, 0)) < 0 ||
              (iJz = SDDS_DefineColumn(&output_defn->SDDSout, "Jz", "J$bz$n", "A/m$a2$n", NULL, NULL, SDDS_FLOAT, 0)) < 0 ||
              (iJr = SDDS_DefineColumn(&output_defn->SDDSout, "Jr", "J$br$n", "A/m$a2$n", NULL, NULL, SDDS_FLOAT, 0)) < 0)
            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
          if (output_defn->separate_imposed_field &&
              ((iEzImposed = SDDS_DefineColumn(&output_defn->SDDSout, "EzImposed", "E$bz,imp$n", "V/m", NULL, NULL, SDDS_FLOAT, 0)) < 0 ||
               (iErImposed = SDDS_DefineColumn(&output_defn->SDDSout, "ErImposed", "E$br,imp$n", "V/m", NULL, NULL, SDDS_FLOAT, 0)) < 0))
            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
        }
      if (EM_problem->modes & FL_TE_FIELDS)
        {
          if ((iEphi = SDDS_DefineColumn(&output_defn->SDDSout, "Ephi", "E$bF$n", "V/m", NULL, NULL, SDDS_FLOAT, 0)) < 0 ||
              (iBr = SDDS_DefineColumn(&output_defn->SDDSout, "Br", "B$br$n", "T", NULL, NULL, SDDS_FLOAT, 0)) < 0 ||
              (iBz = SDDS_DefineColumn(&output_defn->SDDSout, "Bz", "B$bz$n", "T", NULL, NULL, SDDS_FLOAT, 0)) < 0 ||
              (iJphi = SDDS_DefineColumn(&output_defn->SDDSout, "Jphi", "J$bF$n", "A/m$a2$n", NULL, NULL, SDDS_FLOAT, 0)) < 0)
            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
        }
      if (!SDDS_WriteLayout(&output_defn->SDDSout))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
    }
  if (EM_problem->time_step >= output_defn->start_step &&
      (output_defn->step_interval == 0 || (EM_problem->time_step - output_defn->start_step) % output_defn->step_interval == 0))
    {
      output_defn->output_index += 1;
      iRho = SDDS_GetColumnIndex(&output_defn->SDDSout, "rho");
      iEz = SDDS_GetColumnIndex(&output_defn->SDDSout, "Ez");
      iEr = SDDS_GetColumnIndex(&output_defn->SDDSout, "Er");
      iEzImposed = SDDS_GetColumnIndex(&output_defn->SDDSout, "EzImposed");
      iErImposed = SDDS_GetColumnIndex(&output_defn->SDDSout, "ErImposed");
      iJz = SDDS_GetColumnIndex(&output_defn->SDDSout, "Jz");
      iJr = SDDS_GetColumnIndex(&output_defn->SDDSout, "Jr");
      iBphi = SDDS_GetColumnIndex(&output_defn->SDDSout, "Bphi");
      iEphi = SDDS_GetColumnIndex(&output_defn->SDDSout, "Ephi");
      iJphi = SDDS_GetColumnIndex(&output_defn->SDDSout, "Jphi");
      iBr = SDDS_GetColumnIndex(&output_defn->SDDSout, "Br");
      iBz = SDDS_GetColumnIndex(&output_defn->SDDSout, "Bz");
      rIndex = SDDS_GetColumnIndex(&output_defn->SDDSout, "r");
      zIndex = SDDS_GetColumnIndex(&output_defn->SDDSout, "z");
      if (SDDS_NumberOfErrors())
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
      if (!SDDS_StartTable(&output_defn->SDDSout, EM_problem->nz * EM_problem->nr) ||
          !SDDS_SetParameters(&output_defn->SDDSout, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                              "Time", EM_problem->time,
                              "zMinimum", EM_problem->zmin,
                              "zInterval", EM_problem->dz * output_defn->z_interval,
                              "zDimension", (long)ceil(EM_problem->nz / (output_defn->z_interval * 1.0)),
                              "rMinimum", 0.0,
                              "rInterval", EM_problem->dr * output_defn->r_interval,
                              "rDimension", (long)ceil(EM_problem->nr / (output_defn->r_interval * 1.0)),
                              "Timestep", EM_problem->time_step, NULL))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
      row = 0;
      /* output fields at z=iz*dz+zmin and r=ir*dr */
      for (iz = 0; iz < EM_problem->nz; iz += output_defn->z_interval)
        {
          zValue = iz * EM_problem->dz + EM_problem->zmin;
          for (ir = 0; ir < EM_problem->nr; ir += output_defn->r_interval)
            {
              rValue = ir * EM_problem->dr;
              if (!SDDS_SetRowValues(&output_defn->SDDSout, SDDS_BY_INDEX | SDDS_PASS_BY_VALUE, row,
                                     rIndex, rValue, zIndex, zValue,
                                     iRho, EM_problem->Q[iz][ir] / (EM_problem->dz * EM_problem->area_c[ir]), -1))
                {
                  fprintf(stderr, "Problem setting values r, z, rho\n");
                  SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                }
              if (EM_problem->modes & FL_TM_FIELDS)
                {
                  EzImposedValue = (EM_problem->EzImposed[iz][ir] + EM_problem->EzImposed[iz + 1][ir]) / 2;
                  EzValue = (EM_problem->Ez[iz][ir] + EM_problem->Ez[iz + 1][ir]) / 2 +
                    (output_defn->exclude_imposed_field ? 0 : EzImposedValue);
                  ErImposedValue = (EM_problem->ErImposed[iz][ir] + EM_problem->ErImposed[iz + 1][ir]) / 2;
                  ErValue = (EM_problem->Er[iz][ir] + EM_problem->Er[iz][ir + 1]) / 2 +
                    (output_defn->exclude_imposed_field ? 0 : ErImposedValue);
                  JzValue = (EM_problem->Jz[iz][ir] + EM_problem->Jz[iz + 1][ir]) / 2;
                  JrValue = (EM_problem->Jr[iz][ir] + EM_problem->Jr[iz][ir + 1]) / 2;
                  BphiValue = (EM_problem->Bphi2[iz][ir] + EM_problem->Bphi2[iz][ir + 1] +
                               EM_problem->Bphi2[iz + 1][ir] + EM_problem->Bphi2[iz + 1][ir + 1]) /
                    4;
                  if (!SDDS_SetRowValues(&output_defn->SDDSout, SDDS_BY_INDEX | SDDS_PASS_BY_VALUE, row,
                                         iEz, EzValue, iEr, ErValue, iBphi, BphiValue,
                                         iJz, JzValue, iJr, JrValue, -1))
                    {
                      fprintf(stderr, "Problem setting values Ez, Er, Bphi, Jz, Jr\n");
                      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                    }
                  if (output_defn->separate_imposed_field &&
                      !SDDS_SetRowValues(&output_defn->SDDSout, SDDS_BY_INDEX | SDDS_PASS_BY_VALUE, row,
                                         iEzImposed, EzImposedValue, iErImposed, ErImposedValue, -1))
                    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                }
              if (EM_problem->modes & FL_TE_FIELDS)
                {
                  EphiValue = EM_problem->Ephi[iz][ir];
                  JphiValue = EM_problem->Jphi[iz][ir];
                  BrValue = (EM_problem->Br2[iz][ir] + EM_problem->Br2[iz + 1][ir]) / 2;
                  BzValue = (EM_problem->Bz2[iz][ir] + EM_problem->Bz2[iz][ir + 1]) / 2;
                  if (!SDDS_SetRowValues(&output_defn->SDDSout, SDDS_BY_INDEX | SDDS_PASS_BY_VALUE, row,
                                         iEphi, EphiValue, iBr, BrValue, iBz, BzValue,
                                         iJphi, JphiValue, -1))
                    {
                      fprintf(stderr, "Problem setting values Ephi, Br, Bz, Jphi\n");
                      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                    }
                }
              row++;
            }
        }
      if (!SDDS_WriteTable(&output_defn->SDDSout))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
    }
}
