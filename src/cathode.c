/*************************************************************************\
 * Copyright (c) 2002 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2002 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: cathode.c
 * purpose: setup and execution of cathode
 *
 * Michael Borland, 1992 
 */
#include "spiffe.h"
#include "table.h"
#include <time.h>

#define DEBUG 0

void process_cathode_definition(
                                FIELDS *EM_problem,
                                CATHODE **cathodeList,
                                long *cathodes,
                                NAMELIST_TEXT *nl_text)
{
#include "cathode.h"
  CATHODE *cathode;
  long newCathodes, iCathode;
  EMISSION_LOG *emissionLog;

  z_position = inner_radius = outer_radius = current_density = start_time = number_per_step = stop_time =
    discretize_radii = distribution_correction_interval = electrons_per_macroparticle = initial_omega =
    spread_over_dt = determine_temperature = field_emission = 0;
  temperature = work_function = 0;
  stiffness = 1;
  zoned_emission = 1;
  random_number_seed = 987654321;
  profile = profile_time_name = profile_factor_name = NULL;
  field_emission_beta = 1;

  process_namelist(&define_cathode, nl_text);
  print_namelist(stdout, &define_cathode);

  if (z_position < EM_problem->zmin || z_position > EM_problem->zmax)
    bomb("cathode is outside of simulation region", NULL);
  if (inner_radius > outer_radius)
    bomb("cathode has negative radial extent", NULL);
  if (inner_radius == outer_radius && current_density != 0)
    bomb("cathode has zero radial extent--specify electrons_per_macroparticle instead of current_density", NULL);
  if (inner_radius < 0 || outer_radius > EM_problem->rmax)
    bomb("all or part of cathode is outside of simulation region", NULL);
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
  if (temperature > 0 && work_function > 0 && distribution_correction_interval)
    bomb("can't do distribution correction with Richardson-Schottky emission model", NULL);
  if (stop_time && start_time > stop_time)
    bomb("stop_time is non-zero and less than start_time", NULL);
  if (distribution_correction_interval < 0)
    bomb("distribution_correction_interval < 0", NULL);
  if (distribution_correction_interval && discretize_radii)
    bomb("distribution correction and discrete radii are incompatible", NULL);
  if (stiffness <= 0)
    bomb("stiffness <= 0", NULL);
  if (fabs(initial_omega * outer_radius) >= c_mks)
    bomb("initial_omega*outer_radius >= c---unphysical\n", NULL);
  if (add_thermal_velocities && !(temperature > 0 || determine_temperature))
    bomb("add_thermal_velocities is nonzero but the temperature is not positive", NULL);
  if (emission_log)
    emission_log = compose_filename(emission_log, EM_problem->rootname);

  if (random_number_seed == 0)
    {
      random_number_seed = (int)time(NULL);
      random_number_seed = 2 * (random_number_seed / 2) + 1;
      printf("clock-generated random_number_seed = %d\n", random_number_seed);
    }
  random_1(-abs(random_number_seed));

  newCathodes = 1;
  if (field_emission)
    {
      if (inner_radius != 0)
        newCathodes = (outer_radius - inner_radius) / EM_problem->dr;
      else
        newCathodes = (outer_radius - inner_radius + EM_problem->dr / 2) / EM_problem->dr;
    }

  if (!(*cathodeList = SDDS_Realloc(*cathodeList, sizeof(**cathodeList) * (*cathodes + newCathodes))))
    SDDS_Bomb("memory allocation failure");

  emissionLog = NULL;
  for (iCathode = 0; iCathode < newCathodes; iCathode++)
    {
      cathode = *cathodeList + *cathodes + iCathode;
      cathode->z_position = cathode->inner_radius = cathode->outer_radius =
        cathode->current_density = cathode->start_time =
        cathode->start_step = cathode->number_per_step =
        cathode->QoMC = cathode->inner_r2 = cathode->outer_r2 = 0;
      cathode->distribution = cathode->deficit = NULL;
      cathode->area = cathode->rmin = cathode->t_prof = cathode->A_prof = NULL;
      cathode->emissionLog = NULL;

      if (emission_log)
        {
          if (iCathode == 0)
            {
              setUpEmissionLog(&(cathode->emissionLog), emission_log);
              emissionLog = cathode->emissionLog;
            }
          else
            cathode->emissionLog = emissionLog;
        }

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

      cathode->sigmaVMaxwellian = 0;
      if (temperature && add_thermal_velocities)
        cathode->sigmaVMaxwellian = sqrt(temperature * k_boltzmann_mks / me_mks);

      cathode->n_prof_pts = 0;
      if (profile)
        {
          long i;
          SDDS_TABLE SDDSinput;
          if (distribution_correction_interval != 0)
            bomb("can't do distribution correction with an emission profile", NULL);
          if (!SDDS_InitializeInput(&SDDSinput, profile) ||
              SDDS_ReadPage(&SDDSinput) != 1 ||
              !(cathode->t_prof = SDDS_GetColumnInDoubles(&SDDSinput, profile_time_name)) ||
              !(cathode->A_prof = SDDS_GetColumnInDoubles(&SDDSinput, profile_factor_name)))
            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
          if ((cathode->n_prof_pts = SDDS_CountRowsOfInterest(&SDDSinput)) <= 0)
            bomb("no data in profile file", NULL);
          for (i = 1; i < cathode->n_prof_pts; i++)
            if (cathode->t_prof[i - 1] >= cathode->t_prof[i])
              bomb("profile data is not ordered by increasing time values", NULL);
          SDDS_Terminate(&SDDSinput);
        }

      if (!field_emission)
        {
          cathode->inner_radius = inner_radius;
          cathode->outer_radius = outer_radius;
        }
      else
        {
          cathode->outer_radius = (iCathode + 0.5) * EM_problem->dr;
          if ((cathode->inner_radius = cathode->outer_radius - EM_problem->dr) < 0)
            cathode->inner_radius = 0;
          if (iCathode == newCathodes - 1)
            cathode->outer_radius = outer_radius;
          if (cathode->outer_radius > outer_radius)
            cathode->outer_radius = outer_radius;
          if (cathode->inner_radius < inner_radius)
            cathode->inner_radius = inner_radius;
          fprintf(stderr, "Subcathode %ld extends from %e to %e m\n",
                  iCathode, cathode->inner_radius, cathode->outer_radius);
        }
      cathode->z_position = z_position;
      cathode->initial_vz = initial_pz / sqrt(1 + sqr(initial_pz)) * c_mks;
      cathode->initial_omega = initial_omega;
      cathode->current_density = current_density;
      cathode->temperature = temperature;
      cathode->work_function = work_function;
      cathode->field_emission = field_emission;
      cathode->field_emission_beta = field_emission_beta;
      cathode->start_time = start_time;
      cathode->stop_time = stop_time;
      cathode->time_offset = time_offset;
      cathode->autophase = autophase;
      cathode->number_per_step = number_per_step;
      cathode->QoMC = -e_mks / (me_mks * c_mks);
      cathode->correction_interval = distribution_correction_interval;
      cathode->start_step = -1; /* to be set on first call to emit_electrons_from_cathode() */
      cathode->step_emission_started = -1;
      cathode->stiffness = stiffness;
      cathode->electrons_per_macroparticle = electrons_per_macroparticle;
      cathode->spread_over_dt = spread_over_dt;
      cathode->zoned_emission = zoned_emission;
      if ((cathode->discretize_radii = discretize_radii))
        cathode->correction_interval = 0;
      if (cathode->field_emission)
        cathode->zoned_emission = 0;
      if (!cathode->zoned_emission)
        cathode->correction_interval = 0;
      if ((cathode->halton_radix_dt = halton_radix_dt))
        {
          if (!cathode->spread_over_dt)
            bomb("you have nonzero halton_radix_dt but didn't ask to spread_over_dt", NULL);
          if (distribution_correction_interval)
            bomb("you can't use distribution_correction_interval with halton sequence", NULL);
          if (halton_radix_dt < 2)
            bomb("halton_radix_dt<2", NULL);
          cathode->halton_ID_dt = startHaltonSequence(&halton_radix_dt, 0.5);
        }
      if ((cathode->halton_radix_r = halton_radix_r))
        {
          if (distribution_correction_interval)
            bomb("you can't use distribution_correction_interval with halton sequence", NULL);
          if (zoned_emission)
            bomb("you can't use zoned_emission with halton sequence", NULL);
          if (halton_radix_r < 2)
            bomb("halton_radix_r<2", NULL);
          cathode->halton_ID_r = startHaltonSequence(&halton_radix_r, 0.5);
        }
    }

  *cathodes += newCathodes;
}

void emit_electrons_from_cathode(BEAM *beam, CATHODE *cathode, FIELDS *EM_problem, long iCathode)
{
  long ip, np_new;
  long iz, ir, ib, interpCode;
  double total, totalNominal;
  double gamma, dr2, dte, pphi, Ez, Er;
  double area, total_area, r_min, r_max;
  double profile_factor, fnumber;
  double rn;
  double dummy;

  if (cathode->number_per_step == 0 && cathode->electrons_per_macroparticle == 0)
    return;

  iz = (cathode->z_position - EM_problem->zmin) / EM_problem->dz + 1.5;
  ir = (cathode->inner_radius + cathode->outer_radius) / 2 / EM_problem->dr;
  if (iz <= 0)
    iz = 1;
  Ez = EM_problem->constantEz;
  if (EM_problem->modes & FL_TM_FIELDS)
    Ez += EM_problem->Ez[iz][ir] + EM_problem->EzImposed[iz][ir];
  addOffAxisExpansionFields(&Ez, &dummy, &dummy,
                            &dummy, &dummy,
                            cathode->z_position,
                            (cathode->inner_radius + cathode->outer_radius) / 2,
                            EM_problem->time, EM_problem->onAxisField,
                            EM_problem->onAxisFields, 0);
  /*
    if (Ez>=0 && cathode->initial_vz<=0) {
    return;
    }
  */
  if (cathode->start_step == -1)
    {
      if (cathode->current_density <= 1. / sqrt(DBL_MAX))
        return;
      cathode->ref_current_density = cathode->current_density;
      if (cathode->autophase)
        {
          cathode->stop_time = (cathode->stop_time - cathode->start_time) + EM_problem->time + cathode->time_offset;
          cathode->start_time = EM_problem->time + cathode->time_offset;
        }
      if (EM_problem->dti <= 0)
        bomb("integration step size is not defined (emit_electrons_from_cathode)", NULL);
      cathode->inner_r2 = sqr(cathode->inner_radius);
      cathode->outer_r2 = sqr(cathode->outer_radius);
      if ((area = PI * (cathode->outer_r2 - cathode->inner_r2)) != 0)
        {
          if (cathode->electrons_per_macroparticle)
            {
              cathode->number_per_step =
                area * EM_problem->dti * cathode->current_density / (e_mks * cathode->electrons_per_macroparticle);
              if (cathode->number_per_step == 0)
                return;
            }
          else
            cathode->electrons_per_macroparticle =
              area * EM_problem->dti * cathode->current_density / (e_mks * cathode->number_per_step);
        }
      else
        {
          cathode->correction_interval = 0;
          cathode->n_bins = 1;
        }
      cathode->start_step = (cathode->start_time - EM_problem->start_time) / EM_problem->dti;
      if (cathode->stop_time && cathode->stop_time > EM_problem->start_time)
        cathode->stop_step = (cathode->stop_time - EM_problem->start_time) / EM_problem->dti;
      else
        cathode->stop_step = LONG_MAX;
      printf("Each macroparticle will represent %e electrons\n", cathode->electrons_per_macroparticle);
      printf("Initial current density is %e A/m^2\n",
             cathode->current_density);
      cathode->Q = -e_mks * cathode->electrons_per_macroparticle;
      cathode->n_bins = (cathode->outer_radius - cathode->inner_radius) / EM_problem->dr + 1;
      cathode->distribution = tmalloc(sizeof(*cathode->distribution) * cathode->n_bins);
      cathode->deficit = tmalloc(sizeof(*cathode->deficit) * cathode->n_bins);
      cathode->area = tmalloc(sizeof(*cathode->area) * cathode->n_bins);
      cathode->rmin = tmalloc(sizeof(*cathode->rmin) * cathode->n_bins);
      cathode->n_emitted = 0;
      for (ib = cathode->totalArea = 0; ib < cathode->n_bins; ib++)
        {
          ir = cathode->inner_radius / EM_problem->dr + ib;
          if (ir == 0)
            {
              r_max = EM_problem->dr / 2;
              cathode->area[ib] = sqr(r_max) * PI;
              cathode->rmin[ib] = 0;
            }
          else
            {
              r_min = EM_problem->dr * (ir - 0.5);
              r_max = EM_problem->dr * (ir + 0.5);
              if ((cathode->area[ib] = PI * (sqr(r_max) - sqr(r_min))) < 0)
                {
                  fprintf(stderr, "negative area found for cathode annulus, ib=%ld, ir=%ld, dr=%le\n",
                          ib, ir, EM_problem->dr);
                  fprintf(stderr, "inner radius %le -> %le\n", cathode->inner_radius, r_min);
                  fprintf(stderr, "outer radius %le -> %le\n", cathode->outer_radius, r_max);
                  exit(1);
                }
              cathode->rmin[ib] = r_min;
              cathode->totalArea += cathode->area[ib];
            }
#if DEBUG
          if (cathode->area[0] && cathode->zoned_emission)
            printf("cathode annulus %ld area ratio is %e\n", ib, cathode->area[ib] / cathode->area[0]);
#endif
        }
    }

  profile_factor = 1;
  if (cathode->n_prof_pts)
    {
      if (EM_problem->time >= cathode->t_prof[0] && EM_problem->time <= cathode->t_prof[cathode->n_prof_pts - 1])
        profile_factor = interp(cathode->A_prof, cathode->t_prof, cathode->n_prof_pts, EM_problem->time, 1, 1, &interpCode);
      else
        profile_factor = 0;
    }

  if (cathode->start_step > EM_problem->time_step || cathode->stop_step < EM_problem->time_step)
    {
      return;
    }

  iz = (cathode->z_position - EM_problem->zmin) / EM_problem->dz + 1.5;
  ir = cathode->inner_radius / EM_problem->dr;
  if (iz <= 0)
    iz = 1;
  if (cathode->step_emission_started < 0)
    cathode->step_emission_started = EM_problem->time_step;

  dr2 = (cathode->outer_r2 - cathode->inner_r2) / cathode->number_per_step;
  total_area = cathode->totalArea;
  total = totalNominal = cathode->number_per_step *
    (cathode->current_density / cathode->ref_current_density); /* default, without correction */
  np_new = beam->np;
  if (cathode->correction_interval && EM_problem->time_step % cathode->correction_interval == 0)
    {
      total = 0;
      for (ib = 0; ib < cathode->n_bins; ib++)
        {
          ir = cathode->inner_radius / EM_problem->dr + ib;
          if ((cathode->deficit[ib] = ((EM_problem->time_step - cathode->step_emission_started + 1) *
                                       (cathode->area[ib] / total_area) * totalNominal +
                                       0.5) -
               cathode->distribution[ib]) <= 0)
            cathode->deficit[ib] = 0;
          total += cathode->deficit[ib];
#if DEBUG
          printf("distribution[%2ld] = %4ld   deficit[%2ld] = %4ld\n",
                 ib, cathode->distribution[ib], ib, cathode->deficit[ib]);
#endif
        }
    }
  else if (cathode->zoned_emission)
    {
      /* the "deficit" is just the number of particles normally emitted per zone */
      total = 0;
      for (ib = 0; ib < cathode->n_bins; ib++)
        {
          if (total_area)
            {
              fnumber = profile_factor * (cathode->area[ib] / total_area) * totalNominal;
              if (fnumber > 1)
                {
                  cathode->deficit[ib] = fnumber;
                  fnumber -= cathode->deficit[ib];
                }
              else
                cathode->deficit[ib] = 0;
              if (fnumber && fnumber > random_1(1))
                cathode->deficit[ib] += 1;
              total += cathode->deficit[ib];
            }
        }
    }
  else
    {
      fnumber = total * profile_factor;
      total = (long)fnumber;
      fnumber -= (long)fnumber;
      if (fnumber && fnumber > random_1(1))
        total += 1;
    }

  if ((np_new += total) < 0)
    bomb("negative number of particles created!", NULL);

  if (np_new > beam->max_np)
    {
      beam->max_np = 1.1 * (np_new + 1);
      if (beam->max_np > LONG_MAX / 2 || beam->max_np <= 0)
        {
          fprintf(stderr, "Attempted to create %ld particles\n", np_new);
          fprintf(stderr, "number_per_step = %lf\n", cathode->number_per_step);
          fprintf(stderr, "Current density is %e A/m^2\n", cathode->current_density);
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
      if (cathode->correction_interval)
        {
          for (ib = 0; ib < cathode->n_bins; ib++)
            {
              for (ip = beam->np; ip < beam->np + cathode->deficit[ib]; ip++)
                {
                  beam->particleID[ip] = beam->nextParticleID++;
                  beam->rmin[ip] = beam->zmin[ip] = DBL_MAX;
                  beam->rmax[ip] = beam->zmax[ip] = -DBL_MAX;
                  beam->z[ip] = beam->z0[ip] = cathode->z_position;
                  beam->r[ip] = beam->r0[ip] = sqrt(random_1(1) * cathode->area[ib] / PI + sqr(cathode->rmin[ib]));

                  beam->t0[ip] = EM_problem->time;
                  beam->generation[ip] = 0;
                  beam->Q[ip] = cathode->Q;

                  beam->vz[ip] = cathode->initial_vz;
                  beam->vr[ip] = 0;
                  beam->vphi[ip] = cathode->initial_omega * beam->r[ip];
                  if (cathode->sigmaVMaxwellian)
                    addThermalVelocities(beam, ip, cathode->sigmaVMaxwellian);

                  gamma = 1 / sqrt(1 - sqr(beam->vz[ip] / c_mks) - sqr(beam->vr[ip] / c_mks) - sqr(beam->vphi[ip] / c_mks));
                  beam->pz[ip] = beam->vz[ip] / c_mks * gamma;
                  beam->pr[ip] = beam->vr[ip] / c_mks * gamma;
                  beam->r_pphi[ip] = beam->r[ip] * beam->vphi[ip] / c_mks * gamma;

                  beam->status[ip] = PART_ACTIVE;
                  beam->npActive++;

                  if (cathode->spread_over_dt)
                    {
                      iz = (cathode->z_position - EM_problem->zmin) / EM_problem->dz;
                      if (iz <= 0)
                        iz = 1;
                      if (iz > EM_problem->nz)
                        iz = EM_problem->nz;
                      ir = beam->r[ip] / EM_problem->dr;
                      Ez = EM_problem->constantEz;
                      if (EM_problem->modes & FL_TM_FIELDS)
                        Ez = EM_problem->Ez[iz][ir] + EM_problem->EzImposed[iz][ir];
                      ir += 1;
                      iz -= 1;
                      Er = EM_problem->constantEr;
                      if (EM_problem->modes & FL_TM_FIELDS)
                        Er = EM_problem->Er[iz][ir] + EM_problem->ErImposed[iz][ir];
                      dte = EM_problem->dti * (ip - beam->np) / (1.0 * cathode->deficit[ib]);
                      beam->t0[ip] -= dte;
                      pphi = (beam->r[ip] ? beam->r_pphi[ip] / beam->r[ip] : 0);
                      beam->z[ip] += beam->vz[ip] * dte;
                      beam->z0[ip] = beam->z[ip];
                      beam->r[ip] += beam->vr[ip] * dte;
                      beam->r0[ip] = beam->r[ip];
                      beam->pr[ip] += dte * cathode->QoMC * Er / cathode->stiffness;
                      beam->pz[ip] += dte * cathode->QoMC * Ez / cathode->stiffness;
                      gamma = sqrt(sqr(beam->pz[ip]) + sqr(beam->pr[ip]) + sqr(pphi) + 1);
                      beam->vz[ip] = c_mks * beam->pz[ip] / gamma;
                      beam->vr[ip] = c_mks * beam->pr[ip] / gamma;
                    }
                  if (cathode->emissionLog)
                    logCathodeEmission(cathode->emissionLog, beam, ip, iCathode);
                }
              cathode->distribution[ib] += cathode->deficit[ib];
              beam->np = ip;
            }
          if (beam->np != np_new)
            bomb("wrong number of particles emitted for zoned cathode", NULL);
        }
      else if (!cathode->discretize_radii)
        {
          for (ip = beam->np; ip < beam->np + total; ip++)
            {
              beam->particleID[ip] = beam->nextParticleID++;
              beam->rmin[ip] = beam->zmin[ip] = DBL_MAX;
              beam->rmax[ip] = beam->zmax[ip] = -DBL_MAX;
              beam->z[ip] = beam->z0[ip] = cathode->z_position;
              if (cathode->halton_radix_r)
                rn = nextHaltonSequencePoint(cathode->halton_ID_r);
              else
                rn = random_1(1);
              beam->r[ip] = beam->r0[ip] = sqrt(rn *
                                                (sqr(cathode->outer_radius) - sqr(cathode->inner_radius)) +
                                                sqr(cathode->inner_radius));
              beam->t0[ip] = EM_problem->time;
              beam->generation[ip] = 0;
              beam->Q[ip] = cathode->Q;

              beam->vz[ip] = cathode->initial_vz;
              beam->vr[ip] = 0;
              beam->vphi[ip] = cathode->initial_omega * beam->r[ip];
              if (cathode->sigmaVMaxwellian)
                addThermalVelocities(beam, ip, cathode->sigmaVMaxwellian);

              gamma = 1 / sqrt(1 - sqr(beam->vz[ip] / c_mks) - sqr(beam->vr[ip] / c_mks) - sqr(beam->vphi[ip] / c_mks));
              beam->pz[ip] = beam->vz[ip] / c_mks * gamma;
              beam->pr[ip] = beam->vr[ip] / c_mks * gamma;
              beam->r_pphi[ip] = beam->r[ip] * beam->vphi[ip] / c_mks * gamma;

              beam->status[ip] = PART_ACTIVE;
              beam->npActive++;

              if (cathode->spread_over_dt)
                {
                  iz = (cathode->z_position - EM_problem->zmin) / EM_problem->dz;
                  if (iz <= 0)
                    iz = 1;
                  if (iz > EM_problem->nz)
                    iz = EM_problem->nz;
                  ir = beam->r[ip] / EM_problem->dr;
                  Ez = EM_problem->constantEz;
                  if (EM_problem->modes & FL_TM_FIELDS)
                    Ez = EM_problem->Ez[iz][ir] + EM_problem->EzImposed[iz][ir];
                  ir += 1;
                  iz -= 1;
                  Er = EM_problem->constantEr;
                  if (EM_problem->modes & FL_TM_FIELDS)
                    Er = EM_problem->Er[iz][ir] + EM_problem->ErImposed[iz][ir];
                  if (!cathode->halton_radix_dt)
                    dte = EM_problem->dti * (ip - beam->np) / (1.0 * total);
                  else
                    dte = EM_problem->dti * nextHaltonSequencePoint(cathode->halton_ID_dt);
                  beam->t0[ip] -= dte;
                  pphi = (beam->r[ip] ? beam->r_pphi[ip] / beam->r[ip] : 0);
                  beam->z[ip] += beam->vz[ip] * dte;
                  beam->z0[ip] = beam->z[ip];
                  beam->r[ip] += beam->vr[ip] * dte;
                  beam->r0[ip] = beam->r[ip];
                  beam->pr[ip] += dte * cathode->QoMC * Er / cathode->stiffness;
                  beam->pz[ip] += dte * cathode->QoMC * Ez / cathode->stiffness;
                  gamma = sqrt(sqr(beam->pz[ip]) + sqr(beam->pr[ip]) + sqr(pphi) + 1);
                  beam->vz[ip] = c_mks * beam->pz[ip] / gamma;
                  beam->vr[ip] = c_mks * beam->pr[ip] / gamma;
                }
              if (cathode->emissionLog)
                logCathodeEmission(cathode->emissionLog, beam, ip, iCathode);
            }
          beam->np = ip;
          if (beam->np != np_new)
            bomb("wrong number of particles emitted for unzoned cathode", NULL);
        }
      else
        {
          for (ip = beam->np; ip < np_new; ip++)
            {
              beam->particleID[ip] = beam->nextParticleID++;
              beam->rmin[ip] = beam->zmin[ip] = DBL_MAX;
              beam->rmax[ip] = beam->zmax[ip] = -DBL_MAX;
              beam->z[ip] = beam->z0[ip] = cathode->z_position;
              beam->r[ip] = beam->r0[ip] = sqrt(cathode->inner_r2 + (ip - beam->np + 1) * dr2);
              beam->t0[ip] = EM_problem->time;
              beam->generation[ip] = 0;
              beam->Q[ip] = cathode->Q;

              beam->vz[ip] = cathode->initial_vz;
              beam->vr[ip] = 0;
              beam->vphi[ip] = cathode->initial_omega * beam->r[ip];
              if (cathode->sigmaVMaxwellian)
                addThermalVelocities(beam, ip, cathode->sigmaVMaxwellian);

              gamma = 1 / sqrt(1 - sqr(beam->vz[ip] / c_mks) - sqr(beam->vr[ip] / c_mks) - sqr(beam->vphi[ip] / c_mks));
              beam->pz[ip] = beam->vz[ip] / c_mks * gamma;
              beam->pr[ip] = beam->vr[ip] / c_mks * gamma;
              beam->r_pphi[ip] = beam->r[ip] * beam->vphi[ip] / c_mks * gamma;

              beam->status[ip] = PART_ACTIVE;
              beam->npActive++;

              if (cathode->spread_over_dt)
                {
                  iz = (cathode->z_position - EM_problem->zmin) / EM_problem->dz;
                  if (iz <= 0)
                    iz = 1;
                  if (iz > EM_problem->nz)
                    iz = EM_problem->nz;
                  ir = beam->r[ip] / EM_problem->dr;
                  Ez = EM_problem->constantEz;
                  if (EM_problem->modes & FL_TM_FIELDS)
                    Ez = EM_problem->Ez[iz][ir] + EM_problem->EzImposed[iz][ir];
                  ir += 1;
                  iz -= 1;
                  Er = EM_problem->constantEr;
                  if (EM_problem->modes & FL_TM_FIELDS)
                    Er = EM_problem->Er[iz][ir] + EM_problem->ErImposed[iz][ir];
                  dte = EM_problem->dti * (ip - beam->np) / (1.0 * (np_new - beam->np));
                  beam->t0[ip] -= dte;
                  pphi = (beam->r[ip] ? beam->r_pphi[ip] / beam->r[ip] : 0);
                  beam->z[ip] += beam->vz[ip] * dte;
                  beam->z0[ip] = beam->z[ip];
                  beam->r[ip] += beam->vr[ip] * dte;
                  beam->r0[ip] = beam->r[ip];
                  beam->pr[ip] += dte * cathode->QoMC * Er / cathode->stiffness;
                  beam->pz[ip] += dte * cathode->QoMC * Ez / cathode->stiffness;
                  gamma = sqrt(sqr(beam->pz[ip]) + sqr(beam->pr[ip]) + sqr(pphi) + 1);
                  beam->vz[ip] = c_mks * beam->pz[ip] / gamma;
                  beam->vr[ip] = c_mks * beam->pr[ip] / gamma;
                  if (cathode->emissionLog)
                    logCathodeEmission(cathode->emissionLog, beam, ip, iCathode);
                }
            }
        }

      beam->QoMC = cathode->QoMC; /* same for all cathodes (electrons only) */
      beam->np = np_new;
      beam->stiffness = cathode->stiffness;
    }

  if (beam->np > beam->max_np)
    bomb("particle array overflow--programming error (emit_electrons_from_cathode)", NULL);
#if DEBUG
  printf("%ld particles emitted for total of %ld\n", total, beam->np);
#endif
}

void set_current_density(CATHODE *cathode, FIELDS *EM_problem)
{
  double wfCorrection;
  long iz, ir;
  double Ez, dummy;
  iz = (cathode->z_position - EM_problem->zmin) / EM_problem->dz + 1.5;
  ir = (cathode->inner_radius + cathode->outer_radius) / 2 / EM_problem->dr;
  if (iz <= 0)
    iz = 1;
  Ez = EM_problem->constantEz;
  if (EM_problem->modes & FL_TM_FIELDS)
    Ez += EM_problem->Ez[iz][ir] + EM_problem->EzImposed[iz][ir];
  addOffAxisExpansionFields(&Ez, &dummy, &dummy,
                            &dummy, &dummy,
                            cathode->z_position,
                            (cathode->inner_radius + cathode->outer_radius) / 2.0,
                            EM_problem->time, EM_problem->onAxisField,
                            EM_problem->onAxisFields, 0);
  if (Ez >= 0 && cathode->initial_vz <= 0)
    {
      return;
    }
  if (cathode->temperature > 0 && cathode->work_function)
    {
      /* compute current density using Richardson-Schottky equation */
      /* See section 2.4.2.1 of the Handbook of Accelerator Physics and Engineering */
      wfCorrection = 0.012 * sqrt(fabs(Ez * 1e-5));
      if (wfCorrection > cathode->work_function)
        {
          fprintf(stdout, "Warning: Schottky correction exceeds work function.  Set equal to work function.\n");
          wfCorrection = cathode->work_function;
        }
      cathode->current_density =
        60e4 * sqr(cathode->temperature) *
        exp(-(cathode->work_function - wfCorrection) / (cathode->temperature * k_boltzmann_mks / e_mks));
    }
  else if (cathode->field_emission)
    {
      /* compute current density using field emission */
      /* See section 6.12, equation 2 of the Handbook of Accelerator Physics and Engineering */
      cathode->current_density = fieldEmissionCurrentDensity(Ez, cathode->field_emission_beta, cathode->work_function);
    }
}

void logCathodeEmission(EMISSION_LOG *emissionLog, BEAM *beam, long ip, long id)
{
  if (!SDDS_SetRowValues(&(emissionLog->SDDSout), SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE,
                         emissionLog->row,
                         0, beam->t0[ip],
                         1, beam->z0[ip],
                         2, beam->r0[ip],
                         3, beam->vz[ip],
                         4, beam->vr[ip],
                         5, beam->vphi[ip],
                         6, beam->particleID[ip],
                         7, id,
                         -1))
    {
      SDDS_SetError("Problem recording emission log data");
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
    }
  emissionLog->row += 1;
  if (emissionLog->row == emissionLog->maxRows)
    {
      if (!SDDS_WritePage(&(emissionLog->SDDSout)))
        {
          SDDS_SetError("Problem writing emission log data (1)");
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
        }
      emissionLog->row = 0;
      if (!SDDS_StartPage(&(emissionLog->SDDSout), emissionLog->maxRows))
        {
          SDDS_SetError("Problem writing emission log data (2)");
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
        }
    }
}

void finishCathodeEmissionLog(EMISSION_LOG *emissionLog)
{
  if (!SDDS_WritePage(&(emissionLog->SDDSout)) || !SDDS_Terminate(&(emissionLog->SDDSout)))
    {
      SDDS_SetError("Problem writing emission log data");
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
    }
  emissionLog->active = 0;
}

void addThermalVelocities(BEAM *beam, long ip, double sigmaV)
{
  double vx, vy, theta, sin_theta, cos_theta;
  vx = gauss_rn_lim(0.0, sigmaV, 0.0, random_1);
  vy = gauss_rn_lim(0.0, sigmaV, 0.0, random_1);
  theta = random_1(1) * PIx2;
  sin_theta = sin(theta);
  cos_theta = cos(theta);
  beam->vr[ip] += vx * cos_theta + vy * sin_theta;
  beam->vphi[ip] += -vx * sin_theta + vy * cos_theta;
  beam->vz[ip] += fabs(gauss_rn_lim(0.0, sigmaV, 0.0, random_1));
}

void setUpEmissionLog(EMISSION_LOG **emissionLog, char *filename)
{
  SDDS_DATASET *SDDSout;
  *emissionLog = tmalloc(sizeof(**emissionLog));
  SDDSout = &((*emissionLog)->SDDSout);
  (*emissionLog)->row = 0;
  (*emissionLog)->maxRows = 100;
  (*emissionLog)->active = 0;
  if (!SDDS_InitializeOutput(SDDSout, SDDS_BINARY, 0, NULL, NULL, filename) ||
      SDDS_DefineColumn(SDDSout, "t", NULL, "s", NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
      SDDS_DefineColumn(SDDSout, "z", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
      SDDS_DefineColumn(SDDSout, "r", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
      SDDS_DefineColumn(SDDSout, "vz", NULL, "m/s", NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
      SDDS_DefineColumn(SDDSout, "vr", NULL, "m/s", NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
      SDDS_DefineColumn(SDDSout, "vphi", NULL, "m/s", NULL, NULL, SDDS_DOUBLE, 0) < 0 ||
      SDDS_DefineColumn(SDDSout, "particleID", NULL, NULL, NULL, NULL, SDDS_LONG, 0) < 0 ||
      SDDS_DefineColumn(SDDSout, "emitterID", NULL, NULL, NULL, NULL, SDDS_LONG, 0) < 0 ||
      !SDDS_WriteLayout(SDDSout) ||
      !SDDS_StartPage(SDDSout, (*emissionLog)->maxRows))
    {
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
    }
  (*emissionLog)->active = 1;
}

double fieldEmissionCurrentDensity(double E, double beta, double wf)
{
  /* compute current density using field emission */
  /* See section 6.12, equation 2 of the Handbook of Accelerator Physics and Engineering */
  double y, t, v, J;
  E = fabs(E) * beta;
  y = 3.79e-5 * sqrt(E) / wf;
  t = ((-0.0020 * y + 0.0131) * y + 0.1091) * y + 0.9899;
  v = ((0.0445 * y - 0.6782) * y - 0.4091) * y + 1.0437;
  J = 1.54e-6 * sqr(E) / (wf * sqr(t)) * exp(-6.83e9 * v * pow(wf, 1.5) / E);
  return J;
}
