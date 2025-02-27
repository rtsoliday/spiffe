/*************************************************************************\
 * Copyright (c) 2002 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2002 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: test_beam.c
 * purpose: setup and execution of test beam generation
 *
 * Michael Borland, 1992 
 */
#include "spiffe.h"

void process_testbeam_definition(
                                 TEST_BEAM_DEFN *test_beam_defn,
                                 FIELDS *EM_problem,
                                 NAMELIST_TEXT *nl_text)
{
#include "test_beam.h"

  start_time = stop_time = initial_z = initial_r = electrons_per_macroparticle = 0;
  initial_pz = initial_pr = initial_rpphi = 0;
  z_force = r_force = 0;
  z_start_deceleration = r_start_deceleration = DBL_MAX;

  process_namelist(&define_test_beam, nl_text);
  print_namelist(stdout, &define_test_beam);

  if (initial_z < EM_problem->zmin || initial_z > EM_problem->zmax)
    bomb("test beam initial_z is outside of simulation region", NULL);
  if (initial_r < 0 || initial_r > EM_problem->rmax)
    bomb("test beam initial_r is outside of simulation region", NULL);
  if (start_time >= stop_time)
    bomb("start_time >= stop_time for test beam", NULL);
  if (electrons_per_macroparticle <= 0)
    bomb("electrons_per_macroparticle <= 0 for test beam", NULL);
  if (z_force < 0)
    bomb("z_force < 0 for test beam", NULL);

  test_beam_defn->start_time = start_time;
  test_beam_defn->stop_time = stop_time;
  test_beam_defn->initial_z = initial_z;
  test_beam_defn->initial_r = initial_r;
  test_beam_defn->initial_pz = initial_pz;
  test_beam_defn->initial_pr = initial_pr;
  test_beam_defn->initial_rpphi = initial_rpphi;
  test_beam_defn->z_force = z_force;
  test_beam_defn->r_force = r_force;
  test_beam_defn->electrons_per_macroparticle = electrons_per_macroparticle;
  test_beam_defn->start_step = test_beam_defn->stop_step = -1;
  test_beam_defn->z_start_deceleration = z_start_deceleration;
  test_beam_defn->r_start_deceleration = r_start_deceleration;
}

void emit_test_particles(BEAM *beam, TEST_BEAM_DEFN *test_beam_defn, FIELDS *EM_problem)
{
  long ip, np_new;
  double gamma;

  if (test_beam_defn->electrons_per_macroparticle == 0)
    return;

  if (test_beam_defn->start_step == -1)
    {
      if (EM_problem->dti <= 0)
        bomb("integration step size is not defined (emit_electrons_from_cathode)", NULL);
      test_beam_defn->start_step = (test_beam_defn->start_time - EM_problem->start_time) / EM_problem->dti;
      test_beam_defn->stop_step = (test_beam_defn->stop_time - EM_problem->start_time) / EM_problem->dti;
      beam->QoMC = -e_mks / (me_mks * c_mks);
    }

  if (test_beam_defn->start_step <= EM_problem->time_step && test_beam_defn->stop_step >= EM_problem->time_step)
    {
      if ((np_new = beam->np + 1) > beam->max_np)
        {
          beam->max_np = beam->np + 100;
          beam->z = trealloc(beam->z, sizeof(*beam->z) * beam->max_np);
          beam->r = trealloc(beam->r, sizeof(*beam->r) * beam->max_np);
          beam->pz = trealloc(beam->pz, sizeof(*beam->pz) * beam->max_np);
          beam->pr = trealloc(beam->pr, sizeof(*beam->pr) * beam->max_np);
          beam->vz = trealloc(beam->vz, sizeof(*beam->vz) * beam->max_np);
          beam->vr = trealloc(beam->vr, sizeof(*beam->vr) * beam->max_np);
          beam->r0 = trealloc(beam->r0, sizeof(*beam->r0) * beam->max_np);
          beam->t0 = trealloc(beam->t0, sizeof(*beam->t0) * beam->max_np);
          beam->Q = trealloc(beam->Q, sizeof(*beam->Q) * beam->max_np);
          beam->r_pphi = trealloc(beam->r_pphi, sizeof(*beam->r_pphi) * beam->max_np);
          beam->vphi = trealloc(beam->vphi, sizeof(*beam->vphi) * beam->max_np);
        }

      ip = beam->np;
      beam->vz[ip] = beam->vr[ip] = beam->vphi[ip] = 0;
      beam->z[ip] = test_beam_defn->initial_z;
      beam->r[ip] = test_beam_defn->initial_r;
      beam->pz[ip] = test_beam_defn->initial_pz;
      beam->pr[ip] = test_beam_defn->initial_pr;
      beam->r_pphi[ip] = test_beam_defn->initial_rpphi;
      setParticleVelocities(&beam->vz[ip], &beam->vr[ip], &beam->vphi[ip], &gamma,
                            beam->pz[ip], beam->pr[ip], beam->r_pphi[ip], beam->r[ip]);
      beam->Q[ip] = -e_mks * test_beam_defn->electrons_per_macroparticle;
      beam->np = np_new;
    }
}
