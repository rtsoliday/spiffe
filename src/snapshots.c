/*************************************************************************\
 * Copyright (c) 2002 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2002 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: snapshots.c
 * purpose: setup and execution of beam snapshots
 *
 * Michael Borland, 1992 
 */
#include "spiffe.h"

void process_snapshot_definition(
                                 SNAPSHOT *snapshot_defn,
                                 FIELDS *EM_problem,
                                 NAMELIST_TEXT *nl_text)
{
#include "snapshots.h"

  filename = NULL;
  time_interval = start_time = 0;

  process_namelist(&define_snapshots, nl_text);
  print_namelist(stdout, &define_snapshots);

  if (filename == NULL)
    bomb("no filename given for beam snapshots", NULL);
  snapshot_defn->filename = compose_filename(filename, EM_problem->rootname);
  if (time_interval < 0)
    bomb("time_interval < 0 for beam snapshots", NULL);
  snapshot_defn->time_interval = time_interval;
  snapshot_defn->start_time = start_time;

  /* these will be set on the first call to integrate(), when the integration time-step is known */
  snapshot_defn->step_interval = snapshot_defn->start_step = -1;
  snapshot_defn->snapshot_index = -1;
}

#define N_SNAPSHOT_QUANS 9
static char *snapshot_quan[N_SNAPSHOT_QUANS] = {
  "z", "r", "pz", "pr", "pphi", "q", "r0", "z0", "t0"};
static char *snapshot_unit[N_SNAPSHOT_QUANS] = {
  "m", "m", "m$be$nc", "m$be$nc", "m$be$nc", "C", "m", "m", "s"};

void take_beam_snapshot(
                        SNAPSHOT *snapshot_defn,
                        BEAM *beam,
                        FIELDS *EM_problem)
{
  char dump_label[256];
  long ip;
  static double *buffer = NULL;
  static short *status = NULL;

  if (snapshot_defn->filename == NULL)
    return;

  if (snapshot_defn->snapshot_index == -1)
    {
      snapshot_defn->step_interval = snapshot_defn->time_interval / EM_problem->dti + 0.5;
      snapshot_defn->start_step = (snapshot_defn->start_time - EM_problem->start_time) / EM_problem->dti + 0.5;
      snapshot_defn->snapshot_index = 0;
      sprintf(dump_label, "spiffe beam snapshots for run %s", EM_problem->run_name);
      if (!SDDS_InitializeOutput(&snapshot_defn->table, SDDS_BINARY, 0, dump_label, NULL,
                                 snapshot_defn->filename) ||
          !SDDS_DefineSimpleColumns(&snapshot_defn->table, N_SNAPSHOT_QUANS, snapshot_quan, snapshot_unit, SDDS_DOUBLE) ||
          !SDDS_DefineSimpleColumn(&snapshot_defn->table, "generation", "", SDDS_SHORT) ||
          !SDDS_DefineSimpleColumn(&snapshot_defn->table, "particleID", "", SDDS_LONG) ||
          !SDDS_DefineSimpleColumn(&snapshot_defn->table, "status", "", SDDS_SHORT) ||
          SDDS_DefineParameter(&snapshot_defn->table, "Time", NULL, "s", "Time of snapshot", NULL, SDDS_DOUBLE, NULL) < 0 ||
          SDDS_DefineParameter(&snapshot_defn->table, "Timestep", NULL, NULL, "Time-step of snapshot", NULL, SDDS_LONG, NULL) < 0 ||
          !SDDS_WriteLayout(&snapshot_defn->table))
        {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          exit(1);
        }
    }
  if (EM_problem->time_step >= snapshot_defn->start_step &&
      (snapshot_defn->step_interval == 0 || (EM_problem->time_step - snapshot_defn->start_step) % snapshot_defn->step_interval == 0))
    {
      if (!SDDS_StartTable(&snapshot_defn->table, beam->np + 1) ||
          !SDDS_SetParameters(&snapshot_defn->table, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, "Time", EM_problem->time,
                              "Timestep", EM_problem->time_step, NULL))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
      if (!SDDS_SetColumn(&snapshot_defn->table, SDDS_SET_BY_NAME, beam->z, beam->np, "z") ||
          !SDDS_SetColumn(&snapshot_defn->table, SDDS_SET_BY_NAME, beam->r, beam->np, "r") ||
          !SDDS_SetColumn(&snapshot_defn->table, SDDS_SET_BY_NAME, beam->pz, beam->np, "pz") ||
          !SDDS_SetColumn(&snapshot_defn->table, SDDS_SET_BY_NAME, beam->pr, beam->np, "pr") ||
          !SDDS_SetColumn(&snapshot_defn->table, SDDS_SET_BY_NAME, beam->Q, beam->np, "q") ||
          !SDDS_SetColumn(&snapshot_defn->table, SDDS_SET_BY_NAME, beam->r0, beam->np, "r0") ||
          !SDDS_SetColumn(&snapshot_defn->table, SDDS_SET_BY_NAME, beam->z0, beam->np, "z0") ||
          !SDDS_SetColumn(&snapshot_defn->table, SDDS_SET_BY_NAME, beam->t0, beam->np, "t0") ||
          !SDDS_SetColumn(&snapshot_defn->table, SDDS_SET_BY_NAME, beam->generation, beam->np, "generation") ||
          !SDDS_SetColumn(&snapshot_defn->table, SDDS_SET_BY_NAME, beam->particleID, beam->np, "particleID"))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
      buffer = trealloc(buffer, sizeof(*buffer) * beam->np);
      status = trealloc(status, sizeof(*status) * beam->np);
      for (ip = 0; ip < beam->np; ip++)
        {
          buffer[ip] = (beam->r[ip] ? beam->r_pphi[ip] / beam->r[ip] : 0.0);
          status[ip] = beam->status[ip];
        }
      if (!SDDS_SetColumn(&snapshot_defn->table, SDDS_SET_BY_NAME, buffer, beam->np, "pphi") ||
          !SDDS_SetColumn(&snapshot_defn->table, SDDS_SET_BY_NAME, status, beam->np, "status"))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
      if (!SDDS_WriteTable(&snapshot_defn->table))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
      printf("beam snapshot saved to file %s at t=%.16e, step %ld\n",
             snapshot_defn->filename, EM_problem->time, EM_problem->time_step);
      fflush(stdout);
      snapshot_defn->snapshot_index += 1;
    }
}
