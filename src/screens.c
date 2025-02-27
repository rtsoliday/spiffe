/*************************************************************************\
 * Copyright (c) 2002 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2002 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: screens.c
 * purpose: setup beam-imaging screens
 *          screens are executed in advance_particles.c
 *
 * Michael Borland, 1992 
 */
#include "spiffe.h"

static char *screen_quan[N_SCREEN_QUANS] = {
  "r", "pz", "pr", "pphi", "t", "q", "r0", "t0"};
static char *screen_unit[N_SCREEN_QUANS] = {
  "m", "m$be$nc", "m$be$nc", "m$be$nc", "s", "C", "m", "s"};

#define N_SCREEN_DIRECTIONS 2
static char *screen_direction[N_SCREEN_DIRECTIONS] = {
  "forward", "backward"};
static long direction_code[N_SCREEN_DIRECTIONS] = {
  1, -1};

long process_screen_definition(
                               SCREEN **screen_defn,
                               long n_screen_defns,
                               FIELDS *EM_problem,
                               NAMELIST_TEXT *nl_text)
{
#include "screens.h"
  char label[256], name[256];
  long i;

  filename = template = NULL;
  start_time = z_position = delta_z = 0;
  number_of_screens = 1;
  direction = "forward";

  process_namelist(&define_screen, nl_text);
  print_namelist(stdout, &define_screen);

  if (number_of_screens < 1)
    bomb("number_of_screens is invalid", NULL);
  if (filename == NULL && number_of_screens == 1)
    bomb("no filename given for screen", NULL);
  if (template == NULL && number_of_screens > 1)
    bomb("no template given for screens", NULL);
  if (filename != NULL && template != NULL)
    bomb("both filename and template are given--this doesn't make sense", NULL);

  for (i = 0; i < number_of_screens; i++)
    {
      *screen_defn = trealloc(*screen_defn, sizeof(**screen_defn) * (n_screen_defns + 1));
      if (filename)
        {
          (*screen_defn + n_screen_defns)->filename = compose_filename(filename, EM_problem->rootname);
        }
      else
        {
          char *template2;
          template2 = compose_filename(template, EM_problem->rootname);
          sprintf(name, template2, i);
          cp_str(&(*screen_defn + n_screen_defns)->filename, name);
          if (template != template2)
            free(template2);
        }
      (*screen_defn + n_screen_defns)->start_time = start_time;
      (*screen_defn + n_screen_defns)->z_position = z_position;

      if (((*screen_defn + n_screen_defns)->direction = match_string(direction, screen_direction, N_SCREEN_DIRECTIONS, 0)) < 0)
        bomb("unknown screen direction specified", NULL);
      direction = screen_direction[(*screen_defn + n_screen_defns)->direction];
      (*screen_defn + n_screen_defns)->direction = direction_code[(*screen_defn + n_screen_defns)->direction];

      sprintf(label, "spiffe screen output for run %s for particles %s-crossing z=%e m",
              EM_problem->run_name, direction, z_position);
      if (!SDDS_InitializeOutput(&(*screen_defn + n_screen_defns)->table, SDDS_BINARY, 0, label, NULL,
                                 (*screen_defn + n_screen_defns)->filename) ||
          !SDDS_DefineSimpleColumns(&(*screen_defn + n_screen_defns)->table, N_SCREEN_QUANS, screen_quan,
                                    screen_unit, SDDS_DOUBLE) ||
          !SDDS_DefineSimpleColumn(&(*screen_defn + n_screen_defns)->table, "generation", "", SDDS_SHORT) ||
          !SDDS_DefineSimpleColumn(&(*screen_defn + n_screen_defns)->table, "particleID", "", SDDS_LONG) ||
          !SDDS_DefineSimpleColumn(&(*screen_defn + n_screen_defns)->table, "rMin", "m", SDDS_DOUBLE) ||
          !SDDS_DefineSimpleColumn(&(*screen_defn + n_screen_defns)->table, "rMax", "m", SDDS_DOUBLE) ||
          !SDDS_DefineSimpleColumn(&(*screen_defn + n_screen_defns)->table, "zMin", "m", SDDS_DOUBLE) ||
          !SDDS_DefineSimpleColumn(&(*screen_defn + n_screen_defns)->table, "zMax", "m", SDDS_DOUBLE) ||
          SDDS_DefineParameter(&(*screen_defn + n_screen_defns)->table, "direction", NULL, NULL, NULL,
                               NULL, SDDS_STRING, direction) < 0 ||
          SDDS_DefineParameter(&(*screen_defn + n_screen_defns)->table, "Time", NULL, "s", "Time of output",
                               NULL, SDDS_DOUBLE, NULL) < 0 ||
          SDDS_DefineParameter(&(*screen_defn + n_screen_defns)->table, "Timestep", NULL, NULL, "Time-step of output",
                               NULL, SDDS_LONG, NULL) < 0 ||
          !SDDS_WriteLayout(&(*screen_defn + n_screen_defns)->table))
        {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          exit(1);
        }

      (*screen_defn + n_screen_defns)->nb = 0;
      (*screen_defn + n_screen_defns)->max_nb = 0;
      (*screen_defn + n_screen_defns)->buffer = calloc(sizeof(*(*screen_defn + n_screen_defns)->buffer), N_SCREEN_QUANS);
      (*screen_defn + n_screen_defns)->generation = NULL;
      (*screen_defn + n_screen_defns)->particleID = NULL;
      (*screen_defn + n_screen_defns)->rmin =
        (*screen_defn + n_screen_defns)->rmax =
        (*screen_defn + n_screen_defns)->zmin =
        (*screen_defn + n_screen_defns)->zmax = NULL;
      n_screen_defns++;
      z_position += delta_z;
    }
  return (n_screen_defns);
}
