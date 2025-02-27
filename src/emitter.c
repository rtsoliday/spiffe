/*************************************************************************\
 * Copyright (c) 2017 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: emitter.c
 * purpose: setup and execution of emitter regions
 *
 * Michael Borland, 2017
 */
#include "spiffe.h"
#include "table.h"
#include <time.h>

static EMISSION_LOG *emissionLog = NULL;

void process_emitter_definition(
                                FIELDS *EM_problem,
                                EMITTER **emitterList,
                                long *emitters,
                                NAMELIST_TEXT *nl_text)
{
#include "emitter.h"
  EMITTER *emitter;

  current_density = start_time = stop_time = electrons_per_macroparticle =
    determine_temperature = field_emission = 0;
  temperature = work_function = 0;
  random_number_seed = 987654321;
  profile = profile_time_name = profile_factor_name = NULL;
  field_emission_beta = 1;
  halton_radix_r = halton_radix_dt = halton_radix_z = 0;

  process_namelist(&define_emitter, nl_text);
  print_namelist(stdout, &define_emitter);

  if (!field_emission && current_density <= 0 && electrons_per_macroparticle == 0 && temperature <= 0)
    bomb("either field_emission, current_density, electrons_per_macroparticle, or temperature must be non-zero", NULL);
  if (field_emission && (determine_temperature || current_density || profile))
    SDDS_Bomb("field_emission is incompatible with determine_temperature, current_density, and profile");
  if (field_emission && !electrons_per_macroparticle)
    SDDS_Bomb("specify electrons_per_macroparticle with field_emission");
  if (field_emission && !work_function)
    SDDS_Bomb("specify work_function with field_emission");
  if (current_density < 0)
    bomb("the current density should be a positive quantity", NULL);
  if (electrons_per_macroparticle < 0)
    bomb("electrons_per_macroparticle should be a positive quantity", NULL);
  if (temperature < 0)
    bomb("temperature should be a positive quantity", NULL);
  if (stop_time && start_time > stop_time)
    bomb("stop_time is non-zero and less than start_time", NULL);
  if (add_thermal_velocities && !(temperature > 0 || determine_temperature))
    bomb("add_thermal_velocities is nonzero but the temperature is not positive", NULL);
  if (emission_log)
    emission_log = compose_filename(emission_log, EM_problem->rootname);

  if (random_number_seed == 0)
    {
      random_number_seed = (int)time((time_t *)NULL);
      random_number_seed = 2 * (random_number_seed / 2) + 1;
      printf("clock-generated random_number_seed = %d\n", random_number_seed);
    }
  random_2(-abs(random_number_seed));

  *emitterList = SDDS_Realloc(*emitterList, sizeof(**emitterList) * (*emitters + 1));
  emitter = (*emitterList) + (*emitters);
  if ((emitter->nzr = getBoundaryPoints(material_id, &(emitter->iz), &(emitter->ir), &(emitter->direction))) == 0)
    bomb("No surface elements found with the specified material ID", NULL);

  emitter->materialID = material_id;
  emitter->number_per_step = tmalloc(sizeof(*(emitter->number_per_step)) * emitter->nzr);
  emitter->current_density = tmalloc(sizeof(*(emitter->current_density)) * emitter->nzr);

  emitter->initial_v = fabs(initial_p / sqrt(1 + sqr(initial_p)) * c_mks);
  emitter->user_current_density = current_density;
  emitter->temperature = temperature;
  emitter->work_function = work_function;
  emitter->field_emission = field_emission;
  emitter->field_emission_beta = field_emission_beta;
  emitter->start_time = start_time;
  emitter->stop_time = stop_time;
  if ((emitter->electrons_per_macroparticle = electrons_per_macroparticle) == 0)
    bomb("electrons_per_macroparticle must be nonzero", NULL);
  emitter->spread_over_dt = 0;

  emitter->n_emitted = 0;
  emitter->halton_ID_dt = emitter->halton_ID_r = emitter->halton_ID_z = 0;

  emitter->Q = -e_mks * emitter->electrons_per_macroparticle;
  emitter->QoMC = -e_mks / (me_mks * c_mks);
  emitter->start_step = -1; /* to be set on first call to emit_electrons_from_emitter() */
  emitter->step_emission_started = -1;

  emitter->n_prof_pts = 0;
  emitter->t_prof = emitter->A_prof = NULL;
  emitter->offsetFactor = offset_factor;

  emitter->emissionLog = NULL;

  if (determine_temperature)
    {
      long i;
      double lastTemp;
      if (work_function <= 0 || current_density <= 0)
        bomb("to determine temperature, must have positive work function and current density", NULL);
      lastTemp = 1273;
      for (i = 0; i < 100; i++)
        {
          temperature = -e_mks * work_function / k_boltzmann_mks / log(current_density / 60e4 / sqr(lastTemp));
          if (fabs(lastTemp - temperature) < 0.01)
            break;
          lastTemp = temperature;
        }
      fprintf(stdout, "Temperature determined to be %lf K\n", temperature);
    }

  emitter->sigmaVMaxwellian = 0;
  if (temperature && add_thermal_velocities)
    emitter->sigmaVMaxwellian = sqrt(temperature * k_boltzmann_mks / me_mks);

  emitter->n_prof_pts = 0;
  if (profile)
    {
      long i;
      SDDS_TABLE SDDSinput;
      if (!SDDS_InitializeInput(&SDDSinput, profile) ||
          SDDS_ReadPage(&SDDSinput) != 1 ||
          !(emitter->t_prof = SDDS_GetColumnInDoubles(&SDDSinput, profile_time_name)) ||
          !(emitter->A_prof = SDDS_GetColumnInDoubles(&SDDSinput, profile_factor_name)))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
      if ((emitter->n_prof_pts = SDDS_CountRowsOfInterest(&SDDSinput)) <= 0)
        bomb("no data in profile file", NULL);
      for (i = 1; i < emitter->n_prof_pts; i++)
        if (emitter->t_prof[i - 1] >= emitter->t_prof[i])
          bomb("profile data is not ordered by increasing time values", NULL);
      SDDS_Terminate(&SDDSinput);
    }

  if ((emitter->halton_radix_dt = halton_radix_dt))
    {
      if (!emitter->spread_over_dt)
        bomb("you have nonzero halton_radix_dt but didn't ask to spread_over_dt", NULL);
      if (halton_radix_dt < 2)
        bomb("halton_radix_dt<2", NULL);
      emitter->halton_ID_dt = startHaltonSequence(&halton_radix_dt, 0.5);
    }
  if ((emitter->halton_radix_r = halton_radix_r))
    {
      if (halton_radix_r < 2)
        bomb("halton_radix_r<2", NULL);
      emitter->halton_ID_r = startHaltonSequence(&halton_radix_r, 0.5);
    }
  if ((emitter->halton_radix_z = halton_radix_z))
    {
      if (halton_radix_z < 2)
        bomb("halton_radix_z<2", NULL);
      emitter->halton_ID_z = startHaltonSequence(&halton_radix_z, 0.5);
    }

  if (*emitters == 0 && emission_log)
    setUpEmissionLog(&emissionLog, emission_log);

  *emitters += 1;
}

void emit_electrons_from_emitter(BEAM *beam, EMITTER *emitter, FIELDS *EM_problem)
{
  long ip, np_new, ibp;
  long iz, ir, iz1, ir1, interpCode;
  double total, totalNominal;
  double gamma, Ez, Er;
  double r, z, area;
  double profile_factor, fnumber;
  double dummy;

  if (emitter->number_per_step == 0 && emitter->electrons_per_macroparticle == 0)
    return;

  if (EM_problem->dti <= 0)
    bomb("integration step size is not defined (emit_electrons_from_emitter)", NULL);

  if (emitter->start_step > EM_problem->time_step || emitter->stop_step < EM_problem->time_step)
    return;

  profile_factor = 1;
  if (emitter->n_prof_pts)
    {
      if (EM_problem->time >= emitter->t_prof[0] && EM_problem->time <= emitter->t_prof[emitter->n_prof_pts - 1])
        profile_factor = interp(emitter->A_prof, emitter->t_prof, emitter->n_prof_pts, EM_problem->time, 1, 1, &interpCode);
      else
        profile_factor = 0;
    }
  if (profile_factor <= 0)
    return;

  emitter->start_step = (emitter->start_time - EM_problem->start_time) / EM_problem->dti;
  if (emitter->stop_time && emitter->stop_time > EM_problem->start_time)
    emitter->stop_step = (emitter->stop_time - EM_problem->start_time) / EM_problem->dti;
  else
    emitter->stop_step = LONG_MAX;
  if (emitter->start_step > EM_problem->time_step || emitter->stop_step < EM_problem->time_step)
    return;

  for (ibp = 0; ibp < (emitter->nzr - 1); ibp++)
    {
      iz = emitter->iz[ibp];
      ir = emitter->ir[ibp];
      iz1 = emitter->iz[ibp + 1];
      ir1 = emitter->ir[ibp + 1];
      z = emitter->iz[ibp] * EM_problem->dz + EM_problem->zmin;
      r = emitter->ir[ibp] * EM_problem->dr;
      Ez = EM_problem->constantEz;
      Er = EM_problem->constantEr;
      if (EM_problem->modes & FL_TM_FIELDS)
        Ez += EM_problem->Ez[iz + 1][ir] + EM_problem->EzImposed[iz + 1][ir];
      if (EM_problem->modes & FL_TM_FIELDS)
        Er += EM_problem->Er[iz][ir + 1] + EM_problem->ErImposed[iz][ir + 1];
      addOffAxisExpansionFields(&Ez, &Er, &dummy,
                                &dummy, &dummy, z, r,
                                EM_problem->time, EM_problem->onAxisField,
                                EM_problem->onAxisFields, 0);

      emitter->number_per_step[ibp] = 0;
      if (emitter->field_emission)
        emitter->current_density[ibp] =
          fieldEmissionCurrentDensity(sqrt(sqr(Ez) + sqr(Er)), emitter->field_emission_beta, emitter->work_function);
      if (emitter->direction[ibp] == EMITS_RIGHT || emitter->direction[ibp] == EMITS_LEFT)
        area = PI * (sqr(r + EM_problem->dr) - sqr(r));
      else
        area = PIx2 * r * EM_problem->dz;
      if (area > 0 && emitter->electrons_per_macroparticle)
        emitter->number_per_step[ibp] =
          area * EM_problem->dti * emitter->current_density[ibp] / (e_mks * emitter->electrons_per_macroparticle);

      if (emitter->step_emission_started < 0)
        emitter->step_emission_started = EM_problem->time_step;

      total = totalNominal = emitter->number_per_step[ibp];
      np_new = beam->np;

      fnumber = total * profile_factor;
      total = (long)fnumber;
      fnumber -= (long)fnumber;
      if (fnumber && fnumber > random_2(1))
        total += 1;

      if ((np_new += total) < 0)
        bomb("negative number of particles created!", NULL);

      if (np_new > beam->max_np)
        {
          beam->max_np = 1.1 * (np_new + 1);
          if (beam->max_np > LONG_MAX / 2 || beam->max_np <= 0)
            {
              fprintf(stderr, "Attempted to create %ld particles\n", np_new);
              bomb("integer overflow in computing particle array dimension", NULL);
            }
          beam->z = trealloc(beam->z, sizeof(*beam->z) * beam->max_np);
          beam->r = trealloc(beam->r, sizeof(*beam->r) * beam->max_np);
          beam->pz = trealloc(beam->pz, sizeof(*beam->pz) * beam->max_np);
          beam->pr = trealloc(beam->pr, sizeof(*beam->pr) * beam->max_np);
          beam->vz = trealloc(beam->vz, sizeof(*beam->vz) * beam->max_np);
          beam->vr = trealloc(beam->vr, sizeof(*beam->vr) * beam->max_np);
          beam->r0 = trealloc(beam->r0, sizeof(*beam->r0) * beam->max_np);
          beam->z0 = trealloc(beam->z0, sizeof(*beam->z0) * beam->max_np);
          beam->t0 = trealloc(beam->t0, sizeof(*beam->t0) * beam->max_np);
          beam->rmin = trealloc(beam->rmin, sizeof(*beam->rmin) * beam->max_np);
          beam->zmin = trealloc(beam->zmin, sizeof(*beam->zmin) * beam->max_np);
          beam->rmax = trealloc(beam->rmax, sizeof(*beam->rmax) * beam->max_np);
          beam->zmax = trealloc(beam->zmax, sizeof(*beam->zmax) * beam->max_np);
          beam->r_pphi = trealloc(beam->r_pphi, sizeof(*beam->r_pphi) * beam->max_np);
          beam->vphi = trealloc(beam->vphi, sizeof(*beam->vphi) * beam->max_np);
          beam->Q = trealloc(beam->Q, sizeof(*beam->Q) * beam->max_np);
          beam->generation = trealloc(beam->generation, sizeof(*beam->generation) * beam->max_np);
          beam->particleID = trealloc(beam->particleID, sizeof(*beam->particleID) * beam->max_np);
          beam->status = trealloc(beam->status, sizeof(*beam->status) * beam->max_np);
        }

      if (total)
        {
          for (ip = beam->np; ip < beam->np + total; ip++)
            {
              double rn;
              beam->particleID[ip] = beam->nextParticleID++;
              beam->rmin[ip] = beam->zmin[ip] = DBL_MAX;
              beam->rmax[ip] = beam->zmax[ip] = -DBL_MAX;
              beam->z[ip] = z;
              beam->r[ip] = r;
              beam->vz[ip] = beam->vr[ip] = 0;
              if (emitter->direction[ibp] == EMITS_RIGHT || emitter->direction[ibp] == EMITS_LEFT)
                {
                  if (emitter->halton_radix_r)
                    rn = nextHaltonSequencePoint(emitter->halton_ID_r);
                  else
                    rn = random_2(1);
                  if (ir1 > ir)
                    beam->r[ip] = beam->r0[ip] = sqrt(rn * (sqr(r + EM_problem->dr) - sqr(r)) + sqr(r));
                  else
                    beam->r[ip] = beam->r0[ip] = sqrt(rn * (sqr(r + EM_problem->dr) - sqr(r)) + sqr(r - EM_problem->dr));
                  if (emitter->direction[ibp] == EMITS_RIGHT)
                    {
                      beam->vz[ip] = emitter->initial_v;
                      beam->z[ip] += EM_problem->dz * emitter->offsetFactor;
                    }
                  else
                    {
                      beam->vz[ip] = -emitter->initial_v;
                      beam->z[ip] -= EM_problem->dz * emitter->offsetFactor;
                    }
                }
              else
                {
                  if (emitter->halton_radix_z)
                    rn = nextHaltonSequencePoint(emitter->halton_ID_z);
                  else
                    rn = random_2(1);
                  if (iz1 > iz)
                    beam->z[ip] = z + rn * EM_problem->dz;
                  else
                    beam->z[ip] = z - rn * EM_problem->dz;
                  if (emitter->direction[ibp] == EMITS_UP)
                    {
                      beam->vr[ip] = emitter->initial_v;
                      beam->r[ip] += EM_problem->dr * emitter->offsetFactor;
                    }
                  else
                    {
                      beam->vr[ip] = -emitter->initial_v;
                      beam->r[ip] -= EM_problem->dr * emitter->offsetFactor;
                    }
                }

              beam->z0[ip] = beam->z[ip];
              beam->r0[ip] = beam->r[ip];
              beam->t0[ip] = EM_problem->time;
              beam->generation[ip] = 0;
              beam->Q[ip] = emitter->Q;
              beam->vphi[ip] = 0;
              if (emitter->sigmaVMaxwellian)
                addThermalVelocities(beam, ip, emitter->sigmaVMaxwellian);

              gamma = 1 / sqrt(1 - sqr(beam->vz[ip] / c_mks) - sqr(beam->vr[ip] / c_mks) - sqr(beam->vphi[ip] / c_mks));
              beam->pz[ip] = beam->vz[ip] / c_mks * gamma;
              beam->pr[ip] = beam->vr[ip] / c_mks * gamma;
              beam->r_pphi[ip] = beam->r[ip] * beam->vphi[ip] / c_mks * gamma;

              beam->status[ip] = PART_ACTIVE;
              beam->npActive++;

              if (emitter->spread_over_dt)
                {
                }

              if (emissionLog)
                logCathodeEmission(emissionLog, beam, ip, emitter->materialID);
            }
          beam->np = ip;
          if (beam->np != np_new)
            bomb("wrong number of particles emitted for unzoned emitter", NULL);

          beam->QoMC = emitter->QoMC; /* same for all emitters (electrons only) */
          beam->np = np_new;
          beam->stiffness = 1;
        }
    }

  if (beam->np > beam->max_np)
    bomb("particle array overflow--programming error (emit_electrons_from_emitter)", NULL);
#if DEBUG
  printf("%ld particles emitted for total of %ld\n", total, beam->np);
#endif
}

void setEmitterCurrentDensity(EMITTER *emitter, FIELDS *EM_problem)
{
  double wfCorrection;
  long iz, ir;
  double Ez, Er, z, r, dummy;
  long ibp;
  for (ibp = 0; ibp < (emitter->nzr - 1); ibp++)
    {
      emitter->current_density[ibp] = emitter->user_current_density;
      iz = emitter->iz[ibp];
      ir = emitter->ir[ibp];
      z = emitter->iz[ibp] * EM_problem->dz + EM_problem->zmin;
      r = emitter->ir[ibp] * EM_problem->dr;
      Ez = EM_problem->constantEz;
      Er = EM_problem->constantEr;
      if (EM_problem->modes & FL_TM_FIELDS)
        Ez += EM_problem->Ez[iz][ir] + EM_problem->EzImposed[iz][ir];
      if (EM_problem->modes & FL_TM_FIELDS)
        Er += EM_problem->Er[iz][ir] + EM_problem->ErImposed[iz][ir];
      addOffAxisExpansionFields(&Ez, &Er, &dummy,
                                &dummy, &dummy, z, r,
                                EM_problem->time, EM_problem->onAxisField,
                                EM_problem->onAxisFields, 0);

      if (emitter->temperature > 0 && emitter->work_function)
        {
          /* compute current density using Richardson-Schottky equation */
          /* See section 2.4.2.1 of the Handbook of Accelerator Physics and Engineering */
          wfCorrection = 0.012 * sqrt(fabs(Ez * 1e-5));
          if (wfCorrection > emitter->work_function)
            {
              /* fprintf(stdout, "Warning: Schottky correction exceeds work function.  Set equal to work function.\n"); */
              wfCorrection = emitter->work_function;
            }
          emitter->current_density[ibp] =
            60e4 * sqr(emitter->temperature) *
            exp(-(emitter->work_function - wfCorrection) / (emitter->temperature * k_boltzmann_mks / e_mks));
        }
      else if (emitter->field_emission)
        {
          double E;
          E = sqrt(sqr(Ez) + sqr(Er));
          emitter->current_density[ibp] = fieldEmissionCurrentDensity(E, emitter->field_emission_beta, emitter->work_function);
        }
    }
}
