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
 * Michael Borland, 1992-1995
 */
#include "spiffe.h"

void process_field_saving_definition(
                                     FIELD_SAVING *saving_defn,
                                     FIELDS *EM_problem,
                                     NAMELIST_TEXT *nl_text)
{
#include "field_saving.h"

  process_namelist(&define_field_saving, nl_text);
  print_namelist(stdout, &define_field_saving);

  if (!filename)
    bomb("no filename given for field samples", NULL);
  saving_defn->filename = filename;
  if (z_Ez_trigger != -1 && time_interval != -1)
    bomb("Either use time_interval or iz_Ez_trigger, not both", NULL);
  if (time_interval == -1 && z_Ez_trigger == -1)
    bomb("Use one of time_interval or iz_Ez_trigger", NULL);

  saving_defn->time_interval = time_interval;
  saving_defn->start_time = start_time;
  saving_defn->save_before_exiting = save_before_exiting;
  saving_defn->zEzTrigger = z_Ez_trigger;
  saving_defn->EzTriggerLast = 0;

  /* these will be set on the first call to integrate(), when the integration time-step is known */
  saving_defn->step_interval = saving_defn->start_step = -1;
  saving_defn->save_index = -1;
}

void save_fields(
                 FIELD_SAVING *saving_defn,
                 FIELDS *EM_problem,
                 long final_call)
{
  SDDS_TABLE *outTable;
  long iz, ir, iStart, row;
  char *filename;
  long doSave;

  if (saving_defn->filename == NULL)
    return;

  outTable = &saving_defn->outTable;

  if (saving_defn->save_index == -1)
    {
      static char nrBuffer[100], nzBuffer[100];
      static char drBuffer[100], dzBuffer[100];

      sprintf(nrBuffer, "%" PRId32, EM_problem->nr);
      sprintf(nzBuffer, "%" PRId32, EM_problem->nz);
      sprintf(drBuffer, "%21.15e", EM_problem->dr);
      sprintf(dzBuffer, "%21.15e", EM_problem->dz);

      saving_defn->step_interval = saving_defn->time_interval / EM_problem->dti + 0.5;
      saving_defn->start_step = (saving_defn->start_time - EM_problem->start_time) / EM_problem->dti + 0.5;
      saving_defn->save_index = 0;
      filename = compose_filename(saving_defn->filename, EM_problem->rootname);
      if (!SDDS_InitializeOutput(outTable, SDDS_BINARY, 0, NULL, "spiffe field saves",
                                 filename) ||
          SDDS_DefineParameter(outTable, "time", NULL, "s", "time of field save", NULL, SDDS_DOUBLE, NULL) < 0 ||
          SDDS_DefineParameter(outTable, "version", NULL, NULL, "save version", NULL, SDDS_LONG, "4") < 0 ||
          SDDS_DefineParameter1(outTable, "nr", NULL, NULL, "number of r grid lines", NULL, SDDS_LONG,
                                (void *)&EM_problem->nr) < 0 ||
          SDDS_DefineParameter1(outTable, "nz", NULL, NULL, "number of z grid lines", NULL, SDDS_LONG,
                                (void *)&EM_problem->nz) < 0 ||
          SDDS_DefineParameter1(outTable, "dr", NULL, "m", "spacing of r grid lines", NULL, SDDS_DOUBLE,
                                (void *)&EM_problem->dr) < 0 ||
          SDDS_DefineParameter1(outTable, "dz", NULL, "m", "spacing of z grid lines", NULL, SDDS_DOUBLE,
                                (void *)&EM_problem->dz) < 0 ||
          SDDS_DefineParameter1(outTable, "zmin", NULL, "m", NULL, NULL, SDDS_DOUBLE, &EM_problem->zmin) < 0 ||
          SDDS_DefineParameter1(outTable, "zmax", NULL, "m", NULL, NULL, SDDS_DOUBLE, &EM_problem->zmax) < 0 ||
          SDDS_DefineParameter1(outTable, "rmax", NULL, "m", NULL, NULL, SDDS_DOUBLE, &EM_problem->rmax) < 0 ||
          SDDS_DefineParameter1(outTable, "zt", NULL, "m", NULL, NULL, SDDS_DOUBLE, &EM_problem->zt) < 0 ||
          SDDS_DefineParameter(outTable, "TMfields", NULL, NULL, NULL, NULL, SDDS_CHARACTER,
                               EM_problem->modes & FL_TM_FIELDS ? "y" : "n") < 0 ||
          SDDS_DefineParameter(outTable, "TEfields", NULL, NULL, NULL, NULL, SDDS_CHARACTER,
                               EM_problem->modes & FL_TE_FIELDS ? "y" : "n") < 0 ||
          SDDS_DefineColumn(outTable, "metalFlags", NULL, NULL, NULL, NULL, SDDS_LONG, 0) < 0)
        {
          SDDS_SetError("Problem setting up file and parameters for field saving");
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
        }
      if (filename != saving_defn->filename)
        free(filename);
      if (EM_problem->modes & FL_TM_FIELDS &&
          /* the order off these calls should not be changed */
          (SDDS_DefineColumn(outTable, "Ez", NULL, "V/m", NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
           SDDS_DefineColumn(outTable, "Er", NULL, "V/m", NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
           SDDS_DefineColumn(outTable, "Jz", NULL, "A/m$a2$n", NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
           SDDS_DefineColumn(outTable, "Jr", NULL, "A/m$a2$n", NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
           SDDS_DefineColumn(outTable, "Bphi", NULL, "T", NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
           SDDS_DefineColumn(outTable, "Bphi2", NULL, "T", NULL, NULL, SDDS_DOUBLE, 0) < 0))
        {
          SDDS_SetError("Problem defining columns for TM field saving");
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
        }

      if (EM_problem->modes & FL_TE_FIELDS &&
          /* the order off these calls should not be changed */
          (SDDS_DefineColumn(outTable, "Ephi", NULL, "V/m", NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
           SDDS_DefineColumn(outTable, "Jphi", NULL, "A/m$a2$n", NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
           SDDS_DefineColumn(outTable, "Br", NULL, "T", NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
           SDDS_DefineColumn(outTable, "Br2", NULL, "T", NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
           SDDS_DefineColumn(outTable, "Bz", NULL, "T", NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
           SDDS_DefineColumn(outTable, "Bz2", NULL, "T", NULL, NULL, SDDS_DOUBLE, 0) < 0))
        {
          SDDS_SetError("Problem defining columns for TE field saving");
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
        }
      if (!SDDS_WriteLayout(outTable))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
    }
  /* Determine if we should save the fields now */
  doSave = 0;
  iz = (saving_defn->zEzTrigger - EM_problem->zmin) / EM_problem->dz;
  if (final_call && saving_defn->save_before_exiting)
    doSave = 1;
  else if (EM_problem->time_step >= saving_defn->start_step)
    {
      iz = (saving_defn->zEzTrigger - EM_problem->zmin) / EM_problem->dz;
      if (saving_defn->zEzTrigger != -1)
        {
          if (iz >= 0 && iz < EM_problem->nz && SIGN(EM_problem->Ez[iz][0]) != SIGN(saving_defn->EzTriggerLast))
            {
              fprintf(stderr, "Field (z=%le, r=0, iz=%ld) changed from %le to %le---saving\n",
                      saving_defn->zEzTrigger, iz, EM_problem->Ez[iz][0], saving_defn->EzTriggerLast);
              doSave = 1;
            }
        }
      else if (saving_defn->step_interval == 0 || (EM_problem->time_step - saving_defn->start_step) % saving_defn->step_interval == 0)
        doSave = 1;
    }
  if (saving_defn->zEzTrigger != -1 && iz >= 0 && iz < EM_problem->nz)
    saving_defn->EzTriggerLast = EM_problem->Ez[iz][0];

  if (doSave)
    {
      if (!SDDS_StartTable(outTable, (EM_problem->nr + 1) * (EM_problem->nz + 1)) ||
          !SDDS_SetParameters(outTable, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, "time",
                              EM_problem->time, NULL))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);

      if (EM_problem->modes & FL_TM_FIELDS)
        {
          if ((iStart = SDDS_GetColumnIndex(outTable, "Ez")) < 0)
            bomb("problem with field saving--column Ez not present in SDDS layout", NULL);
          for (iz = row = 0; iz <= EM_problem->nz; iz++)
            for (ir = 0; ir <= EM_problem->nr; ir++)
              {
                if (!SDDS_SetRowValues(outTable, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE,
                                       row++,
                                       iStart + 0, EM_problem->Ez[iz][ir],
                                       iStart + 1, EM_problem->Er[iz][ir],
                                       iStart + 2, EM_problem->Jz[iz][ir],
                                       iStart + 3, EM_problem->Jr[iz][ir],
                                       iStart + 4, EM_problem->Bphi[iz][ir],
                                       iStart + 5, EM_problem->Bphi2[iz][ir],
                                       -1))
                  {
                    SDDS_SetError("Problem setting row values for TM field saving");
                    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                  }
              }
        }
      if (EM_problem->modes & FL_TE_FIELDS)
        {
          if ((iStart = SDDS_GetColumnIndex(outTable, "Ephi")) < 0)
            bomb("problem with field saving--column Ephi not present in SDDS layout", NULL);
          for (iz = row = 0; iz <= EM_problem->nz; iz++)
            for (ir = 0; ir <= EM_problem->nr; ir++)
              {
                if (!SDDS_SetRowValues(outTable, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE,
                                       row++,
                                       iStart + 0, EM_problem->Ephi[iz][ir],
                                       iStart + 1, EM_problem->Jphi[iz][ir],
                                       iStart + 2, EM_problem->Br[iz][ir],
                                       iStart + 3, EM_problem->Br2[iz][ir],
                                       iStart + 4, EM_problem->Bz[iz][ir],
                                       iStart + 5, EM_problem->Bz2[iz][ir],
                                       -1))
                  {
                    SDDS_SetError("Problem setting row values for TM field saving");
                    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                  }
              }
        }
      if (!SDDS_WriteTable(outTable))
        {
          SDDS_SetError("Problem writing SDDS table for field saving");
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
        }
      printf("fields saved to file %s at t=%e s\n", saving_defn->filename, EM_problem->time);
      fflush(stdout);
      saving_defn->save_index += 1;
    }
}

int write_array(char **data, long size_of_element, long n1, long n2, FILE *fp)
{
  long i1;

  for (i1 = 0; i1 < n1; i1++)
    if (fwrite(data[i1], size_of_element, n2, fp) != n2)
      return (0);
  return (1);
}
