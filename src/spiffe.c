/*************************************************************************\
 * Copyright (c) 2002 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2002 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: spiffe.c
 *       SPace charge and Integration of Forces For Electrons
 *
 * purpose: main module for particle-in-cell RF gun code
 *
 * Michael Borland, 1992 
 */
#include "spiffe.h"
#include "match_string.h"
#include "rpn.h"
#include <signal.h>
#ifdef SUNOS4
#  include <malloc.h>
#endif

#define DEFINE_GEOMETRY 0
#define DEFINE_ANTENNA 1
#define DEFINE_FIELD_SAMPLING 2
#define DEFINE_FIELD_SAVING 3
#define DEFINE_CATHODE 4
#define LOAD_FIELDS 5
#define INTEGRATE 6
#define TRACE_EXECUTION 7
#define DEFINE_SNAPSHOTS 8
#define DEFINE_SCREEN 9
#define DEFINE_RESISTOR 10
#define DEFINE_TEST_BEAM 11
#define POISSON_CORRECTION 12
#define STOP 13
#define SPACE_CHARGE 14
#define SET_CONSTANT_FIELDS 15
#define SET_UP_TRANSLATION 16
#define LOAD_PARTICLES 17
#define DEFINE_FIELD_OUTPUT 18
#define DEFINE_SOLENOID 19
#define ADD_ON_AXIS_FIELDS 20
#define DEFINE_SECONDARY_EMISSION 21
#define DEFINE_EMITTER 22
#define N_COMMANDS 23
char *command[N_COMMANDS] = {
  "define_geometry",
  "define_antenna",
  "define_field_sampling",
  "define_field_saving",
  "define_cathode",
  "load_fields",
  "integrate",
  "trace",
  "define_snapshots",
  "define_screen",
  "define_resistor",
  "define_test_beam",
  "poisson_correction",
  "stop",
  "space_charge",
  "set_constant_fields",
  "translate",
  "load_particles",
  "define_field_output",
  "define_solenoid",
  "add_on_axis_fields",
  "define_secondary_emission",
  "define_emitter",
};

char *USAGE = "spiffe inputfile\nProgram by M. Borland (This is version 4.10.0, "__DATE__")";

#if defined(CONDOR_COMPILE)
void init_image_with_file_name(char *ckpt_file_name);
#endif

void link_date(void)
{
  fprintf(stdout, "Link date: %s %s, SVN revision: %s\n", __DATE__, __TIME__, SVN_VERSION);
  fflush(stdout);
}
void free_fields_memory(FIELDS *EM_problem);
void free_beam_memory(BEAM *beam);

int main(
         int argc,
         char **argv)
{
  FIELDS EM_problem;        /* the electromagnetic portion of a problem */
  long n_antennae;          /* number of antennae */
  ANTENNA *antenna;         /* pointer to antenna definitions */
  RESISTOR resistor_defns;  /* resitor definitions */
  long n_sample_defns;      /* number of field-sampling definitions */
  SAMPLE_DEFN *sample_defn; /* pointer to field-sampling definitions */
  long n_screen_defns;      /* number of screen definitions */
  SCREEN *screen_defn;      /* pointer to screen definitions */
  SNAPSHOT snapshot_defn;   /* snapshot definition */
  CATHODE *cathode;         /* cathode definition */
  SECONDARY_EMISSION *secondaryEmission = NULL;
  EMITTER *emitter = NULL;
  long nSecondaryDefns = 0, nEmitters = 0;
  long cathodes;
  BEAM beam;                      /* beam properties */
  FIELD_SAVING saving_defn;       /* field-saving definition */
  FIELD_OUTPUT field_output_defn; /* field-output definition (SDDS output) */
  TEST_BEAM_DEFN test_beam_defn;  /* test beam generation parameters */
  long test_beam_defined;
  FILE *fp_in;
  char *inputfile, *ptr;
  char s[1024];
  NAMELIST_TEXT namelist_text;
  long i;

  if (!SDDS_CheckTableStructureSize(sizeof(SDDS_TABLE)))
    {
      fprintf(stderr, "table structure size is inconsistent\n");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }

  if (argc != 2)
    {
      link_date();
      bomb("no input file listed", USAGE);
    }
  inputfile = argv[1];

  if (getenv("RPN_DEFNS"))
    {
      rpn(getenv("RPN_DEFNS"));
      if (rpn_check_error())
        exit(1);
    }

#if defined(CONDOR_COMPILE)
  sprintf(s, "%s.ckpt", inputfile);
  init_image_with_file_name(s);
#endif

  fprintf(stdout, "Running spiffe\n");
  link_date();

#if defined(VAX_VMS) || defined(SUNOS4)
  init_stats();
#endif

  fp_in = fopen_e(inputfile, "r", 0);

  n_antennae = 0;
  antenna = NULL;

  n_sample_defns = 0;
  sample_defn = NULL;
  n_screen_defns = 0;
  screen_defn = NULL;

  saving_defn.filename = NULL;
  saving_defn.time_interval = saving_defn.start_time = 0;
  saving_defn.start_step = 0;

  field_output_defn.filename = NULL;

  beam.np = beam.max_np = beam.npActive = 0;
  beam.z = beam.r = beam.pz = beam.pr = beam.vz = beam.vr = beam.r0 = beam.z0 = beam.t0 =
    beam.r_pphi = beam.vphi = beam.Q = beam.rmin = beam.rmax = beam.zmin = beam.zmax = NULL;
  beam.generation = NULL;
  beam.particleID = NULL;
  beam.nextParticleID = 1;
  beam.status = NULL;
  resistor_defns.n_resistors = 0;

  cathode = NULL;
  cathodes = 0;

  snapshot_defn.filename = NULL;

  SDDS_ZeroMemory(&EM_problem, sizeof(EM_problem));
  EM_problem.nz = EM_problem.nr = EM_problem.time_step = 0;
  EM_problem.dz = EM_problem.dr = EM_problem.zmin = EM_problem.zmax = EM_problem.rmax = 0;
  EM_problem.start_time = EM_problem.dti = EM_problem.time_step = EM_problem.n_time_steps = 0;
  EM_problem.run_name = inputfile;
  SDDS_CopyString(&EM_problem.rootname, inputfile);
  if ((ptr = strchr(EM_problem.rootname, '.')))
    *ptr = 0;

  test_beam_defn.electrons_per_macroparticle = 0; /* indicates that test beam is not defined */

  test_beam_defined = 0;

  while (get_namelist(s, 1024, fp_in))
    {
      scan_namelist(&namelist_text, s);
      switch (match_string(namelist_text.group_name, command, N_COMMANDS, EXACT_MATCH))
        {
        case DEFINE_GEOMETRY:
          process_geometry_definition(&EM_problem, &namelist_text);
          break;
        case DEFINE_ANTENNA:
          process_antenna_definition(&antenna, n_antennae, &EM_problem, &namelist_text);
          n_antennae++;
          break;
        case DEFINE_FIELD_SAMPLING:
          process_sampling_definition(&sample_defn, n_sample_defns, &EM_problem, &namelist_text);
          n_sample_defns++;
          break;
        case DEFINE_FIELD_SAVING:
          process_field_saving_definition(&saving_defn, &EM_problem, &namelist_text);
          break;
        case DEFINE_CATHODE:
          if (test_beam_defined)
            bomb("can't do both cathode and test particles", NULL);
          process_cathode_definition(&EM_problem, &cathode, &cathodes, &namelist_text);
          break;
        case LOAD_FIELDS:
          perform_field_loading(&EM_problem, &namelist_text);
          break;
        case INTEGRATE:
          /* dumpSurfacePointData(EM_problem.rootname, &EM_problem); */
          /* dumpMetalPointData(EM_problem.rootname, &EM_problem); */
          perform_integration(&EM_problem, &beam, cathode, cathodes,
                              &test_beam_defn, antenna, n_antennae, &resistor_defns,
                              sample_defn, n_sample_defns, &saving_defn, &snapshot_defn, screen_defn, n_screen_defns,
                              &field_output_defn, secondaryEmission, nSecondaryDefns, emitter, nEmitters, &namelist_text);
          break;
        case TRACE_EXECUTION:
          fprintf(stderr, "Warning: trace feature is obsolete and no longer available\n");
          break;
        case DEFINE_SNAPSHOTS:
          process_snapshot_definition(&snapshot_defn, &EM_problem, &namelist_text);
          break;
        case DEFINE_SCREEN:
          n_screen_defns = process_screen_definition(&screen_defn, n_screen_defns, &EM_problem, &namelist_text);
          break;
        case DEFINE_RESISTOR:
          process_resistor_definition(&resistor_defns, &EM_problem, &namelist_text);
          break;
        case DEFINE_TEST_BEAM:
          process_testbeam_definition(&test_beam_defn, &EM_problem, &namelist_text);
          if (test_beam_defined)
            bomb("test beam already defined", NULL);
          if (cathodes)
            bomb("can't do both cathode and test beam", NULL);
          break;
        case POISSON_CORRECTION:
          setup_poisson_correction(&EM_problem, &namelist_text);
          break;
        case STOP:
#if defined(VAX_VMS) || defined(SUNOS4)
          report_stats(stdout, "run-time statistics: ");
#endif
          exit(1);
          break;
        case SPACE_CHARGE:
          setup_space_charge(&EM_problem, &namelist_text);
          break;
        case SET_CONSTANT_FIELDS:
          set_constant_fields(&EM_problem, &namelist_text);
          break;
        case SET_UP_TRANSLATION:
          setup_translation(&EM_problem, &namelist_text);
          break;
        case LOAD_PARTICLES:
          perform_particle_loading(&EM_problem, &beam, &namelist_text);
          break;
        case DEFINE_FIELD_OUTPUT:
          process_field_output_definition(&field_output_defn, &EM_problem, &namelist_text);
          break;
        case DEFINE_SOLENOID:
          process_solenoid_definition(&EM_problem, &namelist_text);
          break;
        case ADD_ON_AXIS_FIELDS:
          process_on_axis_field_definition(&EM_problem, &namelist_text);
          break;
        case DEFINE_SECONDARY_EMISSION:
          processSecondaryEmissionDefinition(&secondaryEmission, &nSecondaryDefns, &namelist_text, &EM_problem);
          break;
        case DEFINE_EMITTER:
          process_emitter_definition(&EM_problem, &emitter, &nEmitters, &namelist_text);
          break;
        default:
          printf("unknown namelist %s given.  Known namelists are:\n", namelist_text.group_name);
          for (i = 0; i < N_COMMANDS; i++)
            printf("%s\n", command[i]);
          exit(1);
          break;
        }
      free_namelist_text(&namelist_text);
    }
  if (field_output_defn.filename && !SDDS_Terminate(&field_output_defn.SDDSout))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
  if (snapshot_defn.filename && !SDDS_Terminate(&snapshot_defn.table))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
  for (i = 0; i < n_screen_defns; i++)
    {
      if (screen_defn[i].filename && !SDDS_Terminate(&screen_defn[i].table))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
      if (screen_defn[i].generation)
        free(screen_defn[i].generation);
      if (screen_defn[i].buffer)
        free_zarray_2d((void **)screen_defn[i].buffer, N_SCREEN_QUANS, screen_defn[i].max_nb);
    }
  if (screen_defn)
    free(screen_defn);
  for (i = 0; i < n_sample_defns; i++)
    {
      if (sample_defn[i].filename && !SDDS_Terminate(&sample_defn[i].outTable))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
    }
  if (sample_defn)
    free(sample_defn);
  if (saving_defn.filename && !SDDS_Terminate(&saving_defn.outTable))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
  free_beam_memory(&beam);
  for (i = 0; i < cathodes; i++)
    {
      if (cathode[i].emissionLog)
        finishCathodeEmissionLog(cathode[i].emissionLog);
      if (cathode[i].rmin)
        free(cathode[i].rmin);
      if (cathode[i].area)
        free(cathode[i].area);
      if (cathode[i].deficit)
        free(cathode[i].deficit);
      if (cathode[i].distribution)
        free(cathode[i].distribution);
    }
  if (cathode)
    free(cathode);
  puts("end of input reached");
#if defined(VAX_VMS) || defined(SUNOS4)
  report_stats(stdout, "run-time statistics: ");
#endif
  return (0);
}

char *compose_filename(char *template, char *root_name)
{
  char *ptr;

  if (str_in(template, "%s"))
    {
      ptr = tmalloc(sizeof(char) * (strlen(template) + strlen(root_name) + 1));
      replace_string(ptr, template, "%s", root_name);
      return (ptr);
    }
  else
    return (template);
}

void free_fields_memory(FIELDS *EM_problem)
{
  if (!EM_problem->nr && !EM_problem->nz)
    return;
  if (EM_problem->rootname)
    free(EM_problem->rootname);
  free(EM_problem->area_c);
  free(EM_problem->area_o);
  free_zarray_2d((void **)EM_problem->metal_flags, EM_problem->nz + 1, EM_problem->nr + 1);
  free_zarray_2d((void **)EM_problem->PsiImposed, EM_problem->nz + 1, EM_problem->nr + 1);
  free_zarray_2d((void **)EM_problem->EzImposed, EM_problem->nz + 1, EM_problem->nr + 1);
  free_zarray_2d((void **)EM_problem->ErImposed, EM_problem->nz + 1, EM_problem->nr + 1);
  free_zarray_2d((void **)EM_problem->materialID, EM_problem->nz + 1, EM_problem->nr + 1);
  free_zarray_2d((void **)EM_problem->Q, EM_problem->nz + 1, EM_problem->nr + 1);
  free_zarray_2d((void **)EM_problem->Jz, EM_problem->nz + 1, EM_problem->nr + 1);
  if (EM_problem->modes & FL_TM_FIELDS)
    {
      free_zarray_2d((void **)EM_problem->Ez, EM_problem->nz + 1, EM_problem->nr + 1);
      free_zarray_2d((void **)EM_problem->Er, EM_problem->nz + 1, EM_problem->nr + 1);
      free_zarray_2d((void **)EM_problem->Jr, EM_problem->nz + 1, EM_problem->nr + 1);
      free_zarray_2d((void **)EM_problem->Bphi, EM_problem->nz + 1, EM_problem->nr + 1);
      free_zarray_2d((void **)EM_problem->Bphi2, EM_problem->nz + 1, EM_problem->nr + 1);
      free_zarray_2d((void **)EM_problem->Psi, EM_problem->nz + 1, EM_problem->nr + 1);
    }
  if (EM_problem->modes & FL_TE_FIELDS)
    {
      free_zarray_2d((void **)EM_problem->Ephi, EM_problem->nz + 1, EM_problem->nr + 1);
      free_zarray_2d((void **)EM_problem->Jphi, EM_problem->nz + 1, EM_problem->nr + 1);
      free_zarray_2d((void **)EM_problem->Bz, EM_problem->nz + 1, EM_problem->nr + 1);
      free_zarray_2d((void **)EM_problem->Br, EM_problem->nz + 1, EM_problem->nr + 1);
      free_zarray_2d((void **)EM_problem->Bz2, EM_problem->nz + 1, EM_problem->nr + 1);
      free_zarray_2d((void **)EM_problem->Br2, EM_problem->nz + 1, EM_problem->nr + 1);
    }
}

void free_beam_memory(BEAM *beam)
{
  if (beam->Q)
    free(beam->Q);
  if (beam->z)
    free(beam->z);
  if (beam->r)
    free(beam->r);
  if (beam->pz)
    free(beam->pz);
  if (beam->pr)
    free(beam->pr);
  if (beam->r_pphi)
    free(beam->r_pphi);
  if (beam->vz)
    free(beam->vz);
  if (beam->vr)
    free(beam->vr);
  if (beam->vphi)
    free(beam->vphi);
  if (beam->r0)
    free(beam->r0);
  if (beam->z0)
    free(beam->z0);
  if (beam->t0)
    free(beam->t0);
  if (beam->status)
    free(beam->status);
  if (beam->generation)
    free(beam->generation);
  if (beam->particleID)
    free(beam->particleID);
  if (beam->rmin)
    free(beam->rmin);
  if (beam->rmax)
    free(beam->rmax);
  if (beam->zmin)
    free(beam->zmin);
  if (beam->zmax)
    free(beam->zmax);
}
