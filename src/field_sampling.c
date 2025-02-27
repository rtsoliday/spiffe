/*************************************************************************\
 * Copyright (c) 2002 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2002 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: field_sampling.c
 * purpose: setup and execution of field sampling
 *
 * Michael Borland, 1992 
 */
#include "spiffe.h"

#define SAMPLE_EZ 0
#define SAMPLE_ER 1
#define SAMPLE_EPHI 2
#define SAMPLE_BZ 3
#define SAMPLE_BPHI 4
#define SAMPLE_BR 5
#define SAMPLE_JZ 6
#define SAMPLE_JR 7
#define SAMPLE_PHI 8
#define SAMPLE_Q 9
#define SAMPLE_E 10
#define SAMPLE_B 11
#define N_COMPONENTS 12
static char *component_name[N_COMPONENTS] = {
  "Ez", "Er", "Ephi", "Bz", "Br", "Bphi", "Jz", "Jr", "Phi", "Q", "EAll", "BAll"};
static char *component_symbol[N_COMPONENTS] = {
  "E$bz$n", "E$br$n", "E$b$gf$r$n", "B$bz$n", "B$br$n", "B$b$gf$r$n", "J$bz$n", "J$br$n", "$gF$r", "Q", NULL, NULL};
static char *component_unit[N_COMPONENTS] = {
  "V/m", "V/m", "V/m", "T", "T", "T", "A/m$a2$n", "A/m$a2$n", "V", "C", NULL, NULL};

#define DIRECTION_Z 0
#define DIRECTION_R 1
#define N_DIRECTIONS 2
static char *direction_name[N_DIRECTIONS] = {
  "z", "r"};

void process_sampling_definition(
                                 SAMPLE_DEFN **sample_defn,
                                 long n_sample_defns,
                                 FIELDS *EM_problem,
                                 NAMELIST_TEXT *nl_text)
{
#include "field_sampling.h"
  if (!EM_problem->nr || !EM_problem->nz)
    bomb("you must defined the problem geometry before specifying field sampling", NULL);

  filename = component = direction = NULL;
  min_coord = max_coord = position = time_interval = start_time = time_sequence = 0;

  process_namelist(&define_field_sampling, nl_text);
  print_namelist(stdout, &define_field_sampling);

  *sample_defn = trealloc(*sample_defn, sizeof(**sample_defn) * (n_sample_defns + 1));

  if (filename == NULL)
    bomb("no filename given for field samples", NULL);
  (*sample_defn)[n_sample_defns].filename = compose_filename(filename, EM_problem->rootname);
  if (component == NULL)
    bomb("no field component named for field samples", NULL);
  if (((*sample_defn)[n_sample_defns].component_code = match_string(component, component_name, N_COMPONENTS, 0)) < 0)
    bomb("field component name is not recognized", NULL);
  switch ((*sample_defn)[n_sample_defns].component_code)
    {
    case SAMPLE_EPHI:
    case SAMPLE_BR:
    case SAMPLE_BZ:
      if (!(EM_problem->modes & FL_TE_FIELDS))
        bomb("can't sample TE fields unless they are enabled in define_geometry namelist", NULL);
      break;
    default:
      if (!(EM_problem->modes & FL_TM_FIELDS))
        bomb("can't sample TM fields/potential unless they are enabled in define_geometry namelist", NULL);
      break;
    }
  if (direction == NULL)
    bomb("no direction named for field samples", NULL);
  if (((*sample_defn)[n_sample_defns].direction_code = match_string(direction, direction_name, N_DIRECTIONS, 0)) < 0)
    bomb("direction name is not recognized", NULL);
  if (min_coord > max_coord)
    bomb("min_coord > max_coord for field samples", NULL);
  (*sample_defn)[n_sample_defns].min_coord = min_coord;
  (*sample_defn)[n_sample_defns].max_coord = max_coord;
  (*sample_defn)[n_sample_defns].position = position;
  if (time_interval <= 0)
    bomb("time_interval <=0 for field samples", NULL);
  (*sample_defn)[n_sample_defns].time_interval = time_interval;
  (*sample_defn)[n_sample_defns].start_time = start_time;
  (*sample_defn)[n_sample_defns].time_sequence = time_sequence;

  /* these will be set on the first call to integrate(), when the integration time-step is known */
  (*sample_defn)[n_sample_defns].step_interval = (*sample_defn)[n_sample_defns].start_step = -1;
  (*sample_defn)[n_sample_defns].sample_index = -1;

  (*sample_defn)[n_sample_defns].fileInitialized = 0;
}

void sample_fields(
                   SAMPLE_DEFN *sample_defn,
                   long n_sample_defns,
                   FIELDS *EM_problem)
{
  long i_sample_defn, iz, ir = 0, iz0 = 0, iz1 = 0, ir0, ir1, ifld;
  double value, **data = NULL, **data1;
  char buffer1[256], buffer2[256], buffer3[256];
  double fieldValue[6];
  double fieldSum[6];

  data1 = NULL;
  for (i_sample_defn = 0; i_sample_defn < n_sample_defns; i_sample_defn++)
    {
      if (sample_defn->sample_index == -1)
        {
          sample_defn->step_interval = sample_defn->time_interval / EM_problem->dti + 0.5;
          sample_defn->start_step = (sample_defn->start_time - EM_problem->start_time) / EM_problem->dti + 0.5;
          sample_defn->sample_index = 0;
        }
      switch (sample_defn->component_code)
        {
        case SAMPLE_JZ:
          data = EM_problem->Jz;
          break;
        case SAMPLE_JR:
          data = EM_problem->Jr;
          break;
        case SAMPLE_PHI:
          data = EM_problem->Psi;
          if (EM_problem->PsiImposed)
            data1 = EM_problem->PsiImposed;
          break;
        case SAMPLE_Q:
          data = EM_problem->Q;
          break;
        default:
          break;
        }
      /* This is without a doubt the ugliest, most redundant code I've ever written.
         It does work, however.
      */
      if (!sample_defn->time_sequence)
        {
          if (EM_problem->time_step >= sample_defn->start_step &&
              (sample_defn->step_interval == 0 || (EM_problem->time_step - sample_defn->start_step) % sample_defn->step_interval == 0))
            {
              switch (sample_defn->component_code)
                {
                case SAMPLE_EZ:
                case SAMPLE_BR:
                case SAMPLE_ER:
                case SAMPLE_BZ:
                case SAMPLE_EPHI:
                case SAMPLE_BPHI:
                case SAMPLE_E:
                case SAMPLE_B:
                  switch (sample_defn->direction_code)
                    {
                    case DIRECTION_Z:
                      if ((ir = sample_defn->position / EM_problem->dr) < 0 || ir > EM_problem->nr - 1)
                        bomb("sampling position out of range for Ez, Er, Ephi, Bz, Br, or Bphi samples vs z", NULL);
                      if ((iz0 = (sample_defn->min_coord - EM_problem->zmin) / EM_problem->dz + 0.5) <= 0 || iz0 > EM_problem->nz - 1)
                        iz0 = iz0 <= 0 ? 1 : EM_problem->nz - 1;
                      if ((iz1 = (sample_defn->max_coord - EM_problem->zmin) / EM_problem->dz + 0.5) <= 0 || iz1 > EM_problem->nz - 1)
                        iz1 = iz1 <= 0 ? 1 : EM_problem->nz - 1;
                      if (!sample_defn->fileInitialized)
                        {
                          sprintf(buffer1, "%21.15e", ir * EM_problem->dr);
                          if (!SDDS_InitializeOutput(&sample_defn->outTable, SDDS_BINARY, 0, NULL, NULL,
                                                     sample_defn->filename) ||
                              SDDS_DefineColumn(&sample_defn->outTable, "z", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "rSample", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer1) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "runName", NULL, NULL, NULL, NULL,
                                                   SDDS_STRING, EM_problem->run_name) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "tSample", NULL, "s", NULL, NULL,
                                                   SDDS_DOUBLE, NULL) < 0)
                            {
                              SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                            }
                          if (sample_defn->component_code == SAMPLE_E)
                            {
                              if (SDDS_DefineColumn(&sample_defn->outTable, component_name[0], component_symbol[0], component_unit[0],
                                                    NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
                                  SDDS_DefineColumn(&sample_defn->outTable, component_name[1], component_symbol[1], component_unit[1],
                                                    NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
                                  SDDS_DefineColumn(&sample_defn->outTable, component_name[2], component_symbol[2], component_unit[2],
                                                    NULL, NULL, SDDS_DOUBLE, 0) < 0)
                                SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                            }
                          else if (sample_defn->component_code == SAMPLE_B)
                            {
                              if (SDDS_DefineColumn(&sample_defn->outTable, component_name[3], component_symbol[3], component_unit[3],
                                                    NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
                                  SDDS_DefineColumn(&sample_defn->outTable, component_name[4], component_symbol[4], component_unit[4],
                                                    NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
                                  SDDS_DefineColumn(&sample_defn->outTable, component_name[5], component_symbol[5], component_unit[5],
                                                    NULL, NULL, SDDS_DOUBLE, 0) < 0)
                                SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                            }
                          else if (SDDS_DefineColumn(&sample_defn->outTable,
                                                     component_name[sample_defn->component_code],
                                                     component_symbol[sample_defn->component_code],
                                                     component_unit[sample_defn->component_code],
                                                     NULL, NULL, SDDS_DOUBLE, 0) < 0)
                            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                          if (!SDDS_WriteLayout(&sample_defn->outTable))
                            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                          sample_defn->fileInitialized = 1;
                        }
                      if (!SDDS_StartTable(&sample_defn->outTable, iz1 - iz0 + 1))
                        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                      for (iz = iz0; iz <= iz1; iz++)
                        {
                          interpolateFields(&fieldValue[0], &fieldValue[1], &fieldValue[2],
                                            &fieldValue[3], &fieldValue[4], &fieldValue[5],
                                            EM_problem->zmin + iz * EM_problem->dz,
                                            ir * EM_problem->dr,
                                            EM_problem);
                          if (sample_defn->component_code == SAMPLE_E)
                            {
                              if (!SDDS_SetRowValues(&sample_defn->outTable, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, iz - iz0,
                                                     0, EM_problem->zmin + iz * EM_problem->dz,
                                                     1, fieldValue[0],
                                                     2, fieldValue[1],
                                                     3, fieldValue[2],
                                                     -1))
                                SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                            }
                          else if (sample_defn->component_code == SAMPLE_B)
                            {
                              if (!SDDS_SetRowValues(&sample_defn->outTable, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, iz - iz0,
                                                     0, EM_problem->zmin + iz * EM_problem->dz,
                                                     1, fieldValue[3],
                                                     2, fieldValue[4],
                                                     3, fieldValue[5],
                                                     -1))
                                SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                            }
                          else
                            {
                              if (!SDDS_SetRowValues(&sample_defn->outTable, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, iz - iz0,
                                                     0, EM_problem->zmin + iz * EM_problem->dz,
                                                     1, fieldValue[sample_defn->component_code],
                                                     -1))
                                SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                            }
                        }
                      if (!SDDS_SetParameters(&sample_defn->outTable, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                                              "tSample", EM_problem->time, NULL))
                        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                      if (!SDDS_WriteTable(&sample_defn->outTable))
                        {
                          fprintf(stderr, "Problem writing to %s:\n", sample_defn->filename);
                          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                        }
                      sample_defn->sample_index++;
                      break;
                    case DIRECTION_R:
                      iz = (sample_defn->position - EM_problem->zmin) / EM_problem->dz + 0.5;
                      if (iz <= 0 || iz > EM_problem->nz)
                        bomb("sampling position out of range for Ez, Er, Ephi, Bz, Br, or Bphi samples vs z", NULL);
                      if ((ir0 = sample_defn->min_coord / EM_problem->dr) < 0 || ir0 >= EM_problem->nr)
                        ir0 = ir0 < 0 ? 0 : EM_problem->nr - 1;
                      if ((ir1 = sample_defn->max_coord / EM_problem->dr) < 0 || ir1 >= EM_problem->nr)
                        ir1 = ir1 < 0 ? 0 : EM_problem->nr - 1;
                      if (!sample_defn->fileInitialized)
                        {
                          sprintf(buffer1, "%21.15e", EM_problem->zmin + (iz - 0.5) * EM_problem->dz);
                          if (!SDDS_InitializeOutput(&sample_defn->outTable, SDDS_BINARY, 0, NULL, NULL,
                                                     sample_defn->filename) ||
                              SDDS_DefineColumn(&sample_defn->outTable, "r", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "zSample", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer1) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "runName", NULL, NULL, NULL, NULL,
                                                   SDDS_STRING, EM_problem->run_name) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "tSample", NULL, "s", NULL, NULL,
                                                   SDDS_DOUBLE, NULL) < 0 ||
                              !SDDS_WriteLayout(&sample_defn->outTable))
                            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                          if (sample_defn->component_code == SAMPLE_E)
                            {
                              if (SDDS_DefineColumn(&sample_defn->outTable, component_name[0], component_symbol[0], component_unit[0],
                                                    NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
                                  SDDS_DefineColumn(&sample_defn->outTable, component_name[1], component_symbol[1], component_unit[1],
                                                    NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
                                  SDDS_DefineColumn(&sample_defn->outTable, component_name[2], component_symbol[2], component_unit[2],
                                                    NULL, NULL, SDDS_DOUBLE, 0) < 0)
                                SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                            }
                          else if (sample_defn->component_code == SAMPLE_B)
                            {
                              if (SDDS_DefineColumn(&sample_defn->outTable, component_name[3], component_symbol[3], component_unit[3],
                                                    NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
                                  SDDS_DefineColumn(&sample_defn->outTable, component_name[4], component_symbol[4], component_unit[4],
                                                    NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
                                  SDDS_DefineColumn(&sample_defn->outTable, component_name[5], component_symbol[5], component_unit[5],
                                                    NULL, NULL, SDDS_DOUBLE, 0) < 0)
                                SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                            }
                          else if (SDDS_DefineColumn(&sample_defn->outTable,
                                                     component_name[sample_defn->component_code],
                                                     component_symbol[sample_defn->component_code],
                                                     component_unit[sample_defn->component_code],
                                                     NULL, NULL, SDDS_DOUBLE, 0) < 0)
                            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                          if (!SDDS_WriteLayout(&sample_defn->outTable))
                            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                          sample_defn->fileInitialized = 1;
                        }
                      if (!SDDS_StartTable(&sample_defn->outTable, ir1 - ir0 + 1))
                        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                      for (ir = ir0; ir <= ir1; ir++)
                        {
                          interpolateFields(&fieldValue[0], &fieldValue[1], &fieldValue[2],
                                            &fieldValue[3], &fieldValue[4], &fieldValue[5],
                                            EM_problem->zmin + iz * EM_problem->dz,
                                            ir * EM_problem->dr,
                                            EM_problem);
                          if (sample_defn->component_code == SAMPLE_E)
                            {
                              if (!SDDS_SetRowValues(&sample_defn->outTable, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, iz - iz0,
                                                     0, ir * EM_problem->dr,
                                                     1, fieldValue[0],
                                                     2, fieldValue[1],
                                                     3, fieldValue[2],
                                                     -1))
                                SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                            }
                          else if (sample_defn->component_code == SAMPLE_B)
                            {
                              if (!SDDS_SetRowValues(&sample_defn->outTable, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, iz - iz0,
                                                     0, ir * EM_problem->dr,
                                                     1, fieldValue[3],
                                                     2, fieldValue[4],
                                                     3, fieldValue[5],
                                                     -1))
                                SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                            }
                          else
                            {
                              if (!SDDS_SetRowValues(&sample_defn->outTable, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, iz - iz0,
                                                     0, ir * EM_problem->dr,
                                                     1, fieldValue[sample_defn->component_code],
                                                     -1))
                                SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                            }
                        }
                      if (!SDDS_SetParameters(&sample_defn->outTable, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                                              "tSample", EM_problem->time, NULL))
                        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                      if (!SDDS_WriteTable(&sample_defn->outTable))
                        {
                          fprintf(stderr, "Problem writing to %s:\n", sample_defn->filename);
                          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                        }
                      sample_defn->sample_index++;
                      break;
                    }
                  break;
                case SAMPLE_JZ:
                  switch (sample_defn->direction_code)
                    {
                    case DIRECTION_Z:
                      if ((ir = sample_defn->position / EM_problem->dr) < 0 || ir > EM_problem->nr - 1)
                        bomb("sampling position out of range for Ez/Jz/Br samples vs z", NULL);
                      if ((iz0 = (sample_defn->min_coord - EM_problem->zmin) / EM_problem->dz + 0.5) <= 0 || iz0 > EM_problem->nz - 1)
                        iz0 = iz0 <= 0 ? 1 : EM_problem->nz - 1;
                      if ((iz1 = (sample_defn->max_coord - EM_problem->zmin) / EM_problem->dz + 0.5) <= 0 || iz1 > EM_problem->nz - 1)
                        iz1 = iz1 <= 0 ? 1 : EM_problem->nz - 1;
                      if (!sample_defn->fileInitialized)
                        {
                          sprintf(buffer1, "%21.15e", ir * EM_problem->dr);
                          if (!SDDS_InitializeOutput(&sample_defn->outTable, SDDS_BINARY, 0, NULL, NULL,
                                                     sample_defn->filename) ||
                              SDDS_DefineColumn(&sample_defn->outTable, "z", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
                              SDDS_DefineColumn(&sample_defn->outTable,
                                                component_name[sample_defn->component_code],
                                                component_symbol[sample_defn->component_code],
                                                component_unit[sample_defn->component_code],
                                                NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "rSample", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer1) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "runName", NULL, NULL, NULL, NULL,
                                                   SDDS_STRING, EM_problem->run_name) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "tSample", NULL, "s", NULL, NULL,
                                                   SDDS_DOUBLE, NULL) < 0 ||
                              !SDDS_WriteLayout(&sample_defn->outTable))
                            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                          sample_defn->fileInitialized = 1;
                        }
                      if (!SDDS_StartTable(&sample_defn->outTable, iz1 - iz0 + 1))
                        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                      for (iz = iz0; iz <= iz1; iz++)
                        {
                          if (!SDDS_SetRowValues(&sample_defn->outTable, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, iz - iz0,
                                                 0, EM_problem->zmin + (iz - 0.5) * EM_problem->dz,
                                                 1, data[iz][ir] + (data1 ? data1[iz][ir] : 0.0),
                                                 -1))
                            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                        }
                      if (!SDDS_SetParameters(&sample_defn->outTable, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                                              "tSample", EM_problem->time, NULL))
                        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                      if (!SDDS_WriteTable(&sample_defn->outTable))
                        {
                          fprintf(stderr, "Problem writing to %s:\n", sample_defn->filename);
                          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                        }
                      sample_defn->sample_index++;
                      break;
                    case DIRECTION_R:
                      iz = (sample_defn->position - EM_problem->zmin) / EM_problem->dz + 0.5;
                      if (iz <= 0 || iz > EM_problem->nz)
                        bomb("sampling position out of range for Ez/Jz/Br samples vs r", NULL);
                      if ((ir0 = sample_defn->min_coord / EM_problem->dr) < 0 || ir0 >= EM_problem->nr)
                        ir0 = ir0 < 0 ? 0 : EM_problem->nr - 1;
                      if ((ir1 = sample_defn->max_coord / EM_problem->dr) < 0 || ir1 >= EM_problem->nr)
                        ir1 = ir1 < 0 ? 0 : EM_problem->nr - 1;
                      if (!sample_defn->fileInitialized)
                        {
                          sprintf(buffer1, "%21.15e", EM_problem->zmin + (iz - 0.5) * EM_problem->dz);
                          if (!SDDS_InitializeOutput(&sample_defn->outTable, SDDS_BINARY, 0, NULL, NULL,
                                                     sample_defn->filename) ||
                              SDDS_DefineColumn(&sample_defn->outTable, "r", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
                              SDDS_DefineColumn(&sample_defn->outTable,
                                                component_name[sample_defn->component_code],
                                                component_symbol[sample_defn->component_code],
                                                component_unit[sample_defn->component_code],
                                                NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "zSample", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer1) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "runName", NULL, NULL, NULL, NULL,
                                                   SDDS_STRING, EM_problem->run_name) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "tSample", NULL, "s", NULL, NULL,
                                                   SDDS_DOUBLE, NULL) < 0 ||
                              !SDDS_WriteLayout(&sample_defn->outTable))
                            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                          sample_defn->fileInitialized = 1;
                        }
                      if (!SDDS_StartTable(&sample_defn->outTable, ir1 - ir0 + 1))
                        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                      for (ir = ir0; ir <= ir1; ir++)
                        {
                          if (!SDDS_SetRowValues(&sample_defn->outTable, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, ir - ir0,
                                                 0, EM_problem->dr * ir,
                                                 1, data[iz][ir] + (data1 ? data1[iz][ir] : 0.0),
                                                 -1))
                            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                        }
                      if (!SDDS_SetParameters(&sample_defn->outTable, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                                              "tSample", EM_problem->time, NULL))
                        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                      if (!SDDS_WriteTable(&sample_defn->outTable))
                        {
                          fprintf(stderr, "Problem writing to %s:\n", sample_defn->filename);
                          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                        }
                      sample_defn->sample_index++;
                      break;
                    }
                  break;
                case SAMPLE_JR:
                  switch (sample_defn->direction_code)
                    {
                    case DIRECTION_Z:
                      if ((ir = sample_defn->position / EM_problem->dr + 0.5) < 0 || ir > EM_problem->nr - 1)
                        bomb("sampling position out of range for Er/Jr/Bz samples vs z", NULL);
                      if ((iz0 = (sample_defn->min_coord - EM_problem->zmin) / EM_problem->dz) < 0 || iz0 > EM_problem->nz - 1)
                        iz0 = iz0 < 0 ? 0 : EM_problem->nz - 1;
                      if ((iz1 = (sample_defn->max_coord - EM_problem->zmin) / EM_problem->dz) < 0 || iz1 > EM_problem->nz - 1)
                        iz1 = iz1 < 0 ? 0 : EM_problem->nz - 1;
                      if (!sample_defn->fileInitialized)
                        {
                          sprintf(buffer1, "%21.15e", (ir - 0.5) * EM_problem->dr);
                          if (!SDDS_InitializeOutput(&sample_defn->outTable, SDDS_BINARY, 0, NULL, NULL,
                                                     sample_defn->filename) ||
                              SDDS_DefineColumn(&sample_defn->outTable, "z", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
                              SDDS_DefineColumn(&sample_defn->outTable,
                                                component_name[sample_defn->component_code],
                                                component_symbol[sample_defn->component_code],
                                                component_unit[sample_defn->component_code],
                                                NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "rSample", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer1) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "runName", NULL, NULL, NULL, NULL,
                                                   SDDS_STRING, EM_problem->run_name) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "tSample", NULL, "s", NULL, NULL,
                                                   SDDS_DOUBLE, NULL) < 0 ||
                              !SDDS_WriteLayout(&sample_defn->outTable))
                            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                          sample_defn->fileInitialized = 1;
                        }
                      if (!SDDS_StartTable(&sample_defn->outTable, iz1 - iz0 + 1))
                        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                      for (iz = iz0; iz <= iz1; iz++)
                        {
                          if (!SDDS_SetRowValues(&sample_defn->outTable, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, iz - iz0,
                                                 0, EM_problem->zmin + iz * EM_problem->dz,
                                                 1, data[iz][ir] + (data1 ? data1[iz][ir] : 0.0),
                                                 -1))
                            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                        }
                      if (!SDDS_SetParameters(&sample_defn->outTable, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                                              "tSample", EM_problem->time, NULL))
                        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                      if (!SDDS_WriteTable(&sample_defn->outTable))
                        {
                          fprintf(stderr, "Problem writing to %s:\n", sample_defn->filename);
                          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                        }
                      sample_defn->sample_index++;
                      break;
                    case DIRECTION_R:
                      iz = (sample_defn->position - EM_problem->zmin) / EM_problem->dz;
                      if (iz < 0 || iz >= EM_problem->nz)
                        bomb("sampling position out of range for Er/Jr/Bz samples vs r", NULL);
                      if ((ir0 = sample_defn->min_coord / EM_problem->dr + 0.5) < 1 || ir0 > EM_problem->nr)
                        ir0 = ir0 < 1 ? 1 : EM_problem->nr;
                      if ((ir1 = sample_defn->max_coord / EM_problem->dr + 0.5) < 1 || ir1 > EM_problem->nr)
                        ir1 = ir1 < 1 ? 1 : EM_problem->nr;
                      if (!sample_defn->fileInitialized)
                        {
                          sprintf(buffer1, "%21.15e", EM_problem->zmin + iz * EM_problem->dz);
                          if (!SDDS_InitializeOutput(&sample_defn->outTable, SDDS_BINARY, 0, NULL, NULL,
                                                     sample_defn->filename) ||
                              SDDS_DefineColumn(&sample_defn->outTable, "r", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
                              SDDS_DefineColumn(&sample_defn->outTable,
                                                component_name[sample_defn->component_code],
                                                component_symbol[sample_defn->component_code],
                                                component_unit[sample_defn->component_code],
                                                NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "zSample", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer1) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "runName", NULL, NULL, NULL, NULL,
                                                   SDDS_STRING, EM_problem->run_name) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "tSample", NULL, "s", NULL, NULL,
                                                   SDDS_DOUBLE, NULL) < 0 ||
                              !SDDS_WriteLayout(&sample_defn->outTable))
                            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                          sample_defn->fileInitialized = 1;
                        }
                      if (!SDDS_StartTable(&sample_defn->outTable, ir1 - ir0 + 1))
                        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                      for (ir = ir0; ir <= ir1; ir++)
                        {
                          if (!SDDS_SetRowValues(&sample_defn->outTable, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, ir - ir0,
                                                 0, EM_problem->dr * (ir - 0.5),
                                                 1, data[iz][ir] + (data1 ? data1[iz][ir] : 0.0),
                                                 -1))
                            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                        }
                      if (!SDDS_SetParameters(&sample_defn->outTable, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                                              "tSample", EM_problem->time, NULL))
                        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                      if (!SDDS_WriteTable(&sample_defn->outTable))
                        {
                          fprintf(stderr, "Problem writing to %s:\n", sample_defn->filename);
                          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                        }
                      sample_defn->sample_index++;
                      break;
                    }
                  break;
                case SAMPLE_PHI:
                case SAMPLE_Q:
                  switch (sample_defn->direction_code)
                    {
                    case DIRECTION_Z:
                      if ((ir = sample_defn->position / EM_problem->dr) < 0 || ir > EM_problem->nr - 1)
                        bomb("sampling position out of range for potential/Ephi samples vs z", NULL);
                      if ((iz0 = (sample_defn->min_coord - EM_problem->zmin) / EM_problem->dz) < 0 || iz0 > EM_problem->nz - 1)
                        iz0 = iz0 < 0 ? 0 : EM_problem->nz - 1;
                      if ((iz1 = (sample_defn->max_coord - EM_problem->zmin) / EM_problem->dz) < 0 || iz1 > EM_problem->nz - 1)
                        iz1 = iz1 < 0 ? 0 : EM_problem->nz - 1;
                      if (!sample_defn->fileInitialized)
                        {
                          sprintf(buffer1, "%21.15e", ir * EM_problem->dr);
                          if (!SDDS_InitializeOutput(&sample_defn->outTable, SDDS_BINARY, 0, NULL, NULL,
                                                     sample_defn->filename) ||
                              SDDS_DefineColumn(&sample_defn->outTable, "z", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
                              SDDS_DefineColumn(&sample_defn->outTable,
                                                component_name[sample_defn->component_code],
                                                component_symbol[sample_defn->component_code],
                                                component_unit[sample_defn->component_code],
                                                NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "rSample", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer1) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "runName", NULL, NULL, NULL, NULL,
                                                   SDDS_STRING, EM_problem->run_name) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "tSample", NULL, "s", NULL, NULL,
                                                   SDDS_DOUBLE, NULL) < 0 ||
                              !SDDS_WriteLayout(&sample_defn->outTable))
                            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                          sample_defn->fileInitialized = 1;
                        }
                      if (!SDDS_StartTable(&sample_defn->outTable, iz1 - iz0 + 1))
                        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                      for (iz = iz0; iz <= iz1; iz++)
                        {
                          if (!SDDS_SetRowValues(&sample_defn->outTable, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, iz - iz0,
                                                 0, EM_problem->zmin + iz * EM_problem->dz,
                                                 1, data[iz][ir] + (data1 ? data1[iz][ir] : 0.0),
                                                 -1))
                            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                        }
                      if (!SDDS_SetParameters(&sample_defn->outTable, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                                              "tSample", EM_problem->time, NULL))
                        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                      if (!SDDS_WriteTable(&sample_defn->outTable))
                        {
                          fprintf(stderr, "Problem writing to %s:\n", sample_defn->filename);
                          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                        }
                      sample_defn->sample_index++;
                      break;
                    case DIRECTION_R:
                      iz = (sample_defn->position - EM_problem->zmin) / EM_problem->dz;
                      if (iz < 0 || iz >= EM_problem->nz)
                        bomb("sampling position out of range for potential/Ephi samples vs r", NULL);
                      if ((ir0 = sample_defn->min_coord / EM_problem->dr) < 0 || ir0 >= EM_problem->nr)
                        ir0 = ir0 < 0 ? 0 : EM_problem->nr - 1;
                      if ((ir1 = sample_defn->max_coord / EM_problem->dr) < 0 || ir1 >= EM_problem->nr)
                        ir1 = ir1 < 0 ? 0 : EM_problem->nr - 1;
                      if (!sample_defn->fileInitialized)
                        {
                          sprintf(buffer1, "%21.15e", EM_problem->zmin + iz * EM_problem->dz);
                          if (!SDDS_InitializeOutput(&sample_defn->outTable, SDDS_BINARY, 0, NULL, NULL,
                                                     sample_defn->filename) ||
                              SDDS_DefineColumn(&sample_defn->outTable, "r", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
                              SDDS_DefineColumn(&sample_defn->outTable,
                                                component_name[sample_defn->component_code],
                                                component_symbol[sample_defn->component_code],
                                                component_unit[sample_defn->component_code],
                                                NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "zSample", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer1) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "runName", NULL, NULL, NULL, NULL,
                                                   SDDS_STRING, EM_problem->run_name) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "tSample", NULL, "s", NULL, NULL,
                                                   SDDS_DOUBLE, NULL) < 0 ||
                              !SDDS_WriteLayout(&sample_defn->outTable))
                            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                          sample_defn->fileInitialized = 1;
                        }
                      if (!SDDS_StartTable(&sample_defn->outTable, ir1 - ir0 + 1))
                        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                      for (ir = ir0; ir <= ir1; ir++)
                        {
                          if (!SDDS_SetRowValues(&sample_defn->outTable, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, ir - ir0,
                                                 0, EM_problem->dr * ir,
                                                 1, data[iz][ir] + (data1 ? data1[iz][ir] : 0.0),
                                                 -1))
                            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                        }
                      if (!SDDS_SetParameters(&sample_defn->outTable, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                                              "tSample", EM_problem->time, NULL))
                        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                      if (!SDDS_WriteTable(&sample_defn->outTable))
                        {
                          fprintf(stderr, "Problem writing to %s:\n", sample_defn->filename);
                          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                        }
                      sample_defn->sample_index++;
                      break;
                    }
                  break;
                }
            }
        }
      else
        {
          /* time sequence */
          if (EM_problem->time_step >= sample_defn->start_step &&
              (sample_defn->step_interval == 0 || (EM_problem->time_step - sample_defn->start_step) % sample_defn->step_interval == 0))
            {
              if (!sample_defn->fileInitialized)
                {
                  if (!SDDS_InitializeOutput(&sample_defn->outTable, SDDS_BINARY, 0, NULL, NULL,
                                             sample_defn->filename) ||
                      SDDS_DefineColumn(&sample_defn->outTable, "t", NULL, "s", NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
                      SDDS_DefineParameter(&sample_defn->outTable, "runName", NULL, NULL, NULL, NULL,
                                           SDDS_STRING, EM_problem->run_name) < 0)
                    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                  if (sample_defn->component_code == SAMPLE_E)
                    {
                      if (SDDS_DefineColumn(&sample_defn->outTable, component_name[0], component_symbol[0], component_unit[0],
                                            NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
                          SDDS_DefineColumn(&sample_defn->outTable, component_name[1], component_symbol[1], component_unit[1],
                                            NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
                          SDDS_DefineColumn(&sample_defn->outTable, component_name[2], component_symbol[2], component_unit[2],
                                            NULL, NULL, SDDS_DOUBLE, 0) < 0)
                        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                    }
                  else if (sample_defn->component_code == SAMPLE_B)
                    {
                      if (SDDS_DefineColumn(&sample_defn->outTable, component_name[3], component_symbol[3], component_unit[3],
                                            NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
                          SDDS_DefineColumn(&sample_defn->outTable, component_name[4], component_symbol[4], component_unit[4],
                                            NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
                          SDDS_DefineColumn(&sample_defn->outTable, component_name[5], component_symbol[5], component_unit[5],
                                            NULL, NULL, SDDS_DOUBLE, 0) < 0)
                        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                    }
                  else if (SDDS_DefineColumn(&sample_defn->outTable,
                                             component_name[sample_defn->component_code],
                                             component_symbol[sample_defn->component_code],
                                             component_unit[sample_defn->component_code],
                                             NULL, NULL, SDDS_DOUBLE, 0) < 0)
                    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                }
              switch (sample_defn->component_code)
                {
                case SAMPLE_EZ:
                case SAMPLE_ER:
                case SAMPLE_EPHI:
                case SAMPLE_BZ:
                case SAMPLE_BR:
                case SAMPLE_BPHI:
                case SAMPLE_E:
                case SAMPLE_B:
                  switch (sample_defn->direction_code)
                    {
                    case DIRECTION_Z:
                      if ((ir = sample_defn->position / EM_problem->dr) < 0 || ir > EM_problem->nr - 1)
                        bomb("sampling position out of range for Ez/Jz/Br time-sequence samples vs z", NULL);
                      if ((iz0 = (sample_defn->min_coord - EM_problem->zmin) / EM_problem->dz + 0.5) < 0 || iz0 > EM_problem->nz - 1)
                        bomb("initial sampling coordinate out of range for Ez/Jz/Br samples vs z", NULL);
                      if ((iz1 = (sample_defn->max_coord - EM_problem->zmin) / EM_problem->dz + 0.5) < 0 || iz1 > EM_problem->nz - 1)
                        bomb("final sampling coordinate out of range for Ez/Jz/Br time-sequence samples vs z", NULL);
                      if (iz0 <= 0)
                        iz0 = 1;
                      if (iz1 <= 0)
                        iz1 = 1;
                      if (!sample_defn->fileInitialized)
                        {
                          sprintf(buffer1, "%e", ir * EM_problem->dr);
                          sprintf(buffer2, "%e", iz0 * EM_problem->dz + EM_problem->zmin);
                          sprintf(buffer3, "%e", iz1 * EM_problem->dz + EM_problem->zmin);
                          if (SDDS_DefineParameter(&sample_defn->outTable, "scanDirection", NULL, NULL, NULL, NULL,
                                                   SDDS_STRING, "z") < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "rSample", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer1) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "zStart", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer2) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "zStop", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer3) < 0 ||
                              !SDDS_WriteLayout(&sample_defn->outTable) ||
                              !SDDS_StartTable(&sample_defn->outTable, EM_problem->n_time_steps / sample_defn->step_interval + 1))
                            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                          sample_defn->fileInitialized = 1;
                          sample_defn->sample_index = 0;
                        }
                      for (ifld = 0; ifld < 6; ifld++)
                        fieldSum[ifld] = 0;
                      for (iz = iz0; iz <= iz1; iz++)
                        {
                          interpolateFields(&fieldValue[0], &fieldValue[1], &fieldValue[2],
                                            &fieldValue[3], &fieldValue[4], &fieldValue[5],
                                            EM_problem->zmin + iz * EM_problem->dz,
                                            ir * EM_problem->dr,
                                            EM_problem);
                          for (ifld = 0; ifld < 6; ifld++)
                            fieldSum[ifld] += fieldValue[ifld];
                        }
                      for (ifld = 0; ifld < 6; ifld++)
                        fieldValue[ifld] = fieldSum[ifld] / (iz1 - iz0 + 1);
                      break;
                    case DIRECTION_R:
                      iz = (sample_defn->position - EM_problem->zmin) / EM_problem->dz + 0.5;
                      if (iz <= 0 || iz > EM_problem->nz)
                        bomb("sampling position out of range for Ez/Jz/Br samples vs r", NULL);
                      if ((ir0 = sample_defn->min_coord / EM_problem->dr) < 0 || ir0 >= EM_problem->nr)
                        ir0 = ir0 < 0 ? 0 : EM_problem->nr - 1;
                      if ((ir1 = sample_defn->max_coord / EM_problem->dr) < 0 || ir1 >= EM_problem->nr)
                        ir1 = ir1 < 0 ? 0 : EM_problem->nr - 1;
                      if (!sample_defn->fileInitialized)
                        {
                          sprintf(buffer1, "%e", iz * EM_problem->dz + EM_problem->zmin);
                          sprintf(buffer2, "%e", ir0 * EM_problem->dr);
                          sprintf(buffer3, "%e", ir1 * EM_problem->dr);
                          if (SDDS_DefineParameter(&sample_defn->outTable, "scanDirection", NULL, NULL, NULL, NULL,
                                                   SDDS_STRING, "r") < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "zSample", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer1) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "rStart", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer2) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "rStop", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer3) < 0 ||
                              !SDDS_WriteLayout(&sample_defn->outTable) ||
                              !SDDS_StartTable(&sample_defn->outTable, EM_problem->n_time_steps / sample_defn->step_interval + 1))
                            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                          sample_defn->fileInitialized = 1;
                          sample_defn->sample_index = 0;
                        }
                      for (ifld = 0; ifld < 6; ifld++)
                        fieldSum[ifld] = 0;
                      for (iz = iz0; iz <= iz1; iz++)
                        {
                          interpolateFields(&fieldValue[0], &fieldValue[1], &fieldValue[2],
                                            &fieldValue[3], &fieldValue[4], &fieldValue[5],
                                            EM_problem->zmin + iz * EM_problem->dz,
                                            ir * EM_problem->dr,
                                            EM_problem);
                          for (ifld = 0; ifld < 6; ifld++)
                            fieldSum[ifld] += fieldValue[ifld];
                        }
                      for (ifld = 0; ifld < 6; ifld++)
                        fieldValue[ifld] = fieldSum[ifld] / (iz1 - iz0 + 1);
                      break;
                    }
                  break;
                case SAMPLE_JZ:
                  switch (sample_defn->direction_code)
                    {
                    case DIRECTION_Z:
                      if ((ir = sample_defn->position / EM_problem->dr) < 0 || ir > EM_problem->nr - 1)
                        bomb("sampling position out of range for Ez/Jz/Br time-sequence samples vs z", NULL);
                      if ((iz0 = (sample_defn->min_coord - EM_problem->zmin) / EM_problem->dz + 0.5) < 0 || iz0 > EM_problem->nz - 1)
                        bomb("initial sampling coordinate out of range for Ez/Jz/Br samples vs z", NULL);
                      if ((iz1 = (sample_defn->max_coord - EM_problem->zmin) / EM_problem->dz + 0.5) < 0 || iz1 > EM_problem->nz - 1)
                        bomb("final sampling coordinate out of range for Ez/Jz/Br time-sequence samples vs z", NULL);
                      if (iz0 <= 0)
                        iz0 = 1;
                      if (iz1 <= 0)
                        iz1 = 1;
                      if (!sample_defn->fileInitialized)
                        {
                          sprintf(buffer1, "%e", ir * EM_problem->dr);
                          sprintf(buffer2, "%e", (iz0 - 0.5) * EM_problem->dz + EM_problem->zmin);
                          sprintf(buffer3, "%e", (iz1 - 0.5) * EM_problem->dz + EM_problem->zmin);
                          if (SDDS_DefineParameter(&sample_defn->outTable, "scanDirection", NULL, NULL, NULL, NULL,
                                                   SDDS_STRING, "z") < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "rSample", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer1) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "zStart", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer2) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "zStop", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer3) < 0 ||
                              !SDDS_WriteLayout(&sample_defn->outTable) ||
                              !SDDS_StartTable(&sample_defn->outTable, EM_problem->n_time_steps / sample_defn->step_interval + 1))
                            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                          sample_defn->fileInitialized = 1;
                          sample_defn->sample_index = 0;
                        }
                      value = 0;
                      for (iz = iz0; iz <= iz1; iz++)
                        value += data[iz][ir] + (data1 ? data1[iz][ir] : 0.0);
                      value /= (iz1 - iz0 + 1);
                      break;
                    case DIRECTION_R:
                      iz = (sample_defn->position - EM_problem->zmin) / EM_problem->dz + 0.5;
                      if (iz <= 0 || iz > EM_problem->nz)
                        bomb("sampling position out of range for Ez/Jz/Br samples vs r", NULL);
                      if ((ir0 = sample_defn->min_coord / EM_problem->dr) < 0 || ir0 >= EM_problem->nr)
                        ir0 = ir0 < 0 ? 0 : EM_problem->nr - 1;
                      if ((ir1 = sample_defn->max_coord / EM_problem->dr) < 0 || ir1 >= EM_problem->nr)
                        ir1 = ir1 < 0 ? 0 : EM_problem->nr - 1;
                      if (!sample_defn->fileInitialized)
                        {
                          sprintf(buffer1, "%e", (iz + 0.5) * EM_problem->dz + EM_problem->zmin);
                          sprintf(buffer2, "%e", ir0 * EM_problem->dr);
                          sprintf(buffer3, "%e", ir1 * EM_problem->dr);
                          if (SDDS_DefineParameter(&sample_defn->outTable, "scanDirection", NULL, NULL, NULL, NULL,
                                                   SDDS_STRING, "r") < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "zSample", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer1) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "rStart", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer2) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "rStop", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer3) < 0 ||
                              !SDDS_WriteLayout(&sample_defn->outTable) ||
                              !SDDS_StartTable(&sample_defn->outTable, EM_problem->n_time_steps / sample_defn->step_interval + 1))
                            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                          sample_defn->fileInitialized = 1;
                          sample_defn->sample_index = 0;
                        }
                      value = 0;
                      for (ir = ir0; ir <= ir1; ir++)
                        value += data[iz][ir] + (data1 ? data1[iz][ir] : 0.0);
                      value /= (ir1 - ir0 + 1);
                      break;
                    }
                  break;
                case SAMPLE_JR:
                  switch (sample_defn->direction_code)
                    {
                    case DIRECTION_Z:
                      if ((ir = sample_defn->position / EM_problem->dr + 0.5) < 0 || ir > EM_problem->nr - 1)
                        bomb("sampling position out of range for Er/Jr/Bz samples vs z", NULL);
                      if ((iz0 = (sample_defn->min_coord - EM_problem->zmin) / EM_problem->dz) < 0 || iz0 > EM_problem->nz - 1)
                        bomb("initial sampling coordinate out of range for Er/Jr/Bz samples vs z", NULL);
                      if ((iz1 = (sample_defn->max_coord - EM_problem->zmin) / EM_problem->dz) < 0 || iz1 > EM_problem->nz - 1)
                        bomb("final sampling coordinate out of range for Er/Jr/Bz samples vs z", NULL);
                      if (!sample_defn->fileInitialized)
                        {
                          sprintf(buffer1, "%e", (ir - 0.5) * EM_problem->dr);
                          sprintf(buffer2, "%e", iz0 * EM_problem->dz + EM_problem->zmin);
                          sprintf(buffer3, "%e", iz1 * EM_problem->dz + EM_problem->zmin);
                          if (SDDS_DefineParameter(&sample_defn->outTable, "scanDirection", NULL, NULL, NULL, NULL,
                                                   SDDS_STRING, "z") < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "rSample", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer1) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "zStart", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer2) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "zStop", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer3) < 0 ||
                              !SDDS_WriteLayout(&sample_defn->outTable) ||
                              !SDDS_StartTable(&sample_defn->outTable, EM_problem->n_time_steps / sample_defn->step_interval + 1))
                            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                          sample_defn->fileInitialized = 1;
                          sample_defn->sample_index = 0;
                        }
                      value = 0;
                      for (iz = iz0; iz <= iz1; iz++)
                        value += data[iz][ir] + (data1 ? data1[iz][ir] : 0.0);
                      value /= (iz1 - iz0 + 1);
                      break;
                    case DIRECTION_R:
                      iz = (sample_defn->position - EM_problem->zmin) / EM_problem->dz;
                      if (iz < 0 || iz >= EM_problem->nz)
                        bomb("sampling position out of range for Er/Jr/Bz samples vs r", NULL);
                      if ((ir0 = sample_defn->min_coord / EM_problem->dr + 0.5) < 1 || ir0 > EM_problem->nr)
                        ir0 = ir0 < 1 ? 1 : EM_problem->nr;
                      if ((ir1 = sample_defn->max_coord / EM_problem->dr + 0.5) < 1 || ir1 > EM_problem->nr)
                        ir1 = ir1 < 1 ? 1 : EM_problem->nr;
                      if (!sample_defn->fileInitialized)
                        {
                          sprintf(buffer1, "%e", iz * EM_problem->dz + EM_problem->zmin);
                          sprintf(buffer2, "%e", (ir0 + 0.5) * EM_problem->dr);
                          sprintf(buffer3, "%e", (ir1 + 0.5) * EM_problem->dr);
                          if (SDDS_DefineParameter(&sample_defn->outTable, "scanDirection", NULL, NULL, NULL, NULL,
                                                   SDDS_STRING, "r") < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "zSample", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer1) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "rStart", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer2) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "rStop", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer3) < 0 ||
                              !SDDS_WriteLayout(&sample_defn->outTable) ||
                              !SDDS_StartTable(&sample_defn->outTable, EM_problem->n_time_steps / sample_defn->step_interval + 1))
                            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                          sample_defn->fileInitialized = 1;
                          sample_defn->sample_index = 0;
                        }
                      value = 0;
                      for (ir = ir0; ir <= ir1; ir++)
                        value += data[iz][ir] + (data1 ? data1[iz][ir] : 0.0);
                      value /= (ir1 - ir0 + 1);
                      break;
                    }
                  break;
                case SAMPLE_PHI:
                case SAMPLE_Q:
                  switch (sample_defn->direction_code)
                    {
                    case DIRECTION_Z:
                      if ((ir = sample_defn->position / EM_problem->dr) < 0 || ir > EM_problem->nr - 1)
                        bomb("sampling position out of range for potential/Ephi samples vs z", NULL);
                      if ((iz0 = (sample_defn->min_coord - EM_problem->zmin) / EM_problem->dz) < 0 || iz0 > EM_problem->nz - 1)
                        bomb("initial sampling coordinate out of range for potential/Ephi samples vs z", NULL);
                      if ((iz1 = (sample_defn->max_coord - EM_problem->zmin) / EM_problem->dz) < 0 || iz1 > EM_problem->nz - 1)
                        bomb("final sampling coordinate out of range for potential/Ephi samples vs z", NULL);
                      if (!sample_defn->fileInitialized)
                        {
                          sprintf(buffer1, "%e", ir * EM_problem->dr);
                          sprintf(buffer2, "%e", iz0 * EM_problem->dz + EM_problem->zmin);
                          sprintf(buffer3, "%e", iz1 * EM_problem->dz + EM_problem->zmin);
                          if (SDDS_DefineParameter(&sample_defn->outTable, "scanDirection", NULL, NULL, NULL, NULL,
                                                   SDDS_STRING, "z") < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "rSample", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer1) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "zStart", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer2) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "zStop", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer3) < 0 ||
                              !SDDS_WriteLayout(&sample_defn->outTable) ||
                              !SDDS_StartTable(&sample_defn->outTable, EM_problem->n_time_steps / sample_defn->step_interval + 1))
                            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                          sample_defn->fileInitialized = 1;
                          sample_defn->sample_index = 0;
                        }
                      value = 0;
                      for (iz = iz0; iz <= iz1; iz++)
                        value += data[iz][ir] + (data1 ? data1[iz][ir] : 0.0);
                      value /= (iz1 - iz0 + 1);
                      break;
                    case DIRECTION_R:
                      iz = (sample_defn->position - EM_problem->zmin) / EM_problem->dz;
                      if (iz < 0 || iz >= EM_problem->nz)
                        bomb("sampling position out of range for potential/Ephi samples vs r", NULL);
                      if ((ir0 = sample_defn->min_coord / EM_problem->dr) < 0 || ir0 >= EM_problem->nr)
                        ir0 = ir0 < 0 ? 0 : EM_problem->nr - 1;
                      if ((ir1 = sample_defn->max_coord / EM_problem->dr) < 0 || ir1 >= EM_problem->nr)
                        ir1 = ir1 < 0 ? 0 : EM_problem->nr - 1;
                      if (!sample_defn->fileInitialized)
                        {
                          sprintf(buffer1, "%e", iz * EM_problem->dz + EM_problem->zmin);
                          sprintf(buffer2, "%e", ir0 * EM_problem->dr);
                          sprintf(buffer3, "%e", ir1 * EM_problem->dr);
                          if (SDDS_DefineParameter(&sample_defn->outTable, "scanDirection", NULL, NULL, NULL, NULL,
                                                   SDDS_STRING, "r") < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "zSample", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer1) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "rStart", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer2) < 0 ||
                              SDDS_DefineParameter(&sample_defn->outTable, "rStop", NULL, "m", NULL, NULL,
                                                   SDDS_DOUBLE, buffer3) < 0 ||
                              !SDDS_WriteLayout(&sample_defn->outTable) ||
                              !SDDS_StartTable(&sample_defn->outTable, EM_problem->n_time_steps / sample_defn->step_interval + 1))
                            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                          sample_defn->fileInitialized = 1;
                          sample_defn->sample_index = 0;
                        }
                      value = 0;
                      for (ir = ir0; ir <= ir1; ir++)
                        value += data[iz][ir] + (data1 ? data1[iz][ir] : 0.0);
                      value /= (ir1 - ir0 + 1);
                      break;
                    }
                }
              switch (sample_defn->component_code)
                {
                case SAMPLE_EZ:
                case SAMPLE_ER:
                case SAMPLE_EPHI:
                case SAMPLE_BZ:
                case SAMPLE_BR:
                case SAMPLE_BPHI:
                  if (!SDDS_SetRowValues(&sample_defn->outTable, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE,
                                         sample_defn->sample_index,
                                         0, EM_problem->time_step * EM_problem->dti + EM_problem->start_time,
                                         1, fieldValue[sample_defn->component_code],
                                         -1) ||
                      !SDDS_UpdateTable(&sample_defn->outTable))
                    {
                      fprintf(stderr, "Problem setting/writing SDDS data for time sequence field samples, file %s.\n",
                              sample_defn->filename);
                      fprintf(stderr, "Row is %ld\n", sample_defn->sample_index);
                      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                    }
                  break;
                case SAMPLE_E:
                  if (!SDDS_SetRowValues(&sample_defn->outTable, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE,
                                         sample_defn->sample_index,
                                         0, EM_problem->time_step * EM_problem->dti + EM_problem->start_time,
                                         1, fieldValue[0],
                                         2, fieldValue[1],
                                         3, fieldValue[2],
                                         -1) ||
                      !SDDS_UpdateTable(&sample_defn->outTable))
                    {
                      fprintf(stderr, "Problem setting/writing SDDS data for time sequence field samples, file %s.\n",
                              sample_defn->filename);
                      fprintf(stderr, "Row is %ld\n", sample_defn->sample_index);
                      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                    }
                  break;
                case SAMPLE_B:
                  break;
                default:
                  if (!SDDS_SetRowValues(&sample_defn->outTable, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE,
                                         sample_defn->sample_index,
                                         0, EM_problem->time_step * EM_problem->dti + EM_problem->start_time,
                                         1, fieldValue[3],
                                         2, fieldValue[4],
                                         3, fieldValue[5],
                                         -1) ||
                      !SDDS_UpdateTable(&sample_defn->outTable))
                    {
                      fprintf(stderr, "Problem setting/writing SDDS data for time sequence field samples, file %s.\n",
                              sample_defn->filename);
                      fprintf(stderr, "Row is %ld\n", sample_defn->sample_index);
                      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
                    }
                }
              sample_defn->sample_index++;
            }
        }
      sample_defn++;
    }
}
