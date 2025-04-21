/*************************************************************************\
 * Copyright (c) 2002 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2002 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: integrate.c
 * purpose: do integration
 *
 * Michael Borland, 1992 
 */
#include "spiffe.h"

void calculate_areas(FIELDS *EM_problem);

long isqr(long i)
{
  return (i * i);
}

void perform_integration(
                         FIELDS *EM_problem,
                         BEAM *beam,
                         CATHODE *cathode,
                         long cathodes,
                         TEST_BEAM_DEFN *test_beam_defn,
                         ANTENNA *antenna,
                         long n_antennae,
                         RESISTOR *resistor,
                         SAMPLE_DEFN *sample_defn,
                         long n_sample_defns,
                         FIELD_SAVING *field_saving_defn,
                         SNAPSHOT *snapshot_defn,
                         SCREEN *screen_defn,
                         long n_screen_defns,
                         FIELD_OUTPUT *field_output_defn,
                         SECONDARY_EMISSION *secondaryEmission,
                         long nSecondaryDefns,
                         EMITTER *emitter,
                         long nEmitters,
                         NAMELIST_TEXT *nl_text /* unparsed integration parameters in namelist text format */
                         )
{
#include "integrate.h"
  double dtC, Jm1, Jm2;
  long is_reload, beam_seen;
  SDDS_DATASET SDDSstatus, SDDSlost;
  long outputRow = 0, icat;

  dt_integration = start_time = finish_time = space_charge = 0;
  check_divergence = 0;
  status_interval = -1;
  smoothing_parameter = 0;
  J_filter_multiplier = 0;
  terminate_on_total_loss = SG_smoothing_halfwidth = 0;
  lost_particles = status_output = NULL;
  z_forces_start = -1;

  process_namelist(&integrate, nl_text);
  print_namelist(stdout, &integrate);

  if (dt_integration <= 0)
    bomb("dt_integration <= 0", NULL);
  if (finish_time < start_time)
    bomb("finish_time < start_time", NULL);
  if (smoothing_parameter < 0 || smoothing_parameter > 1)
    bomb("smoothing_parameter must be >= 0 and < 1", NULL);
  if (J_filter_multiplier < 0 || J_filter_multiplier >= 1)
    bomb("J_filter_multiplier must be >=0 and < 1", NULL);
  if (EM_problem->dti)
    is_reload = 1;
  else
    is_reload = 0;
  if (space_charge)
    EM_problem->modes |= FL_SPACE_CHARGE;
  if (status_output)
    setUpStatusOutputFile(&SDDSstatus, status_output, EM_problem->rootname);
  if (lost_particles)
    setUpLostParticleFile(&SDDSlost, lost_particles, EM_problem->rootname);

  EM_problem->imposedEFieldRampTime = imposed_Efield_ramp_time;
  EM_problem->imposedEFieldFlatTopTime = imposed_Efield_flat_top_time;
  EM_problem->dti = dt_integration;
  EM_problem->status_interval = status_interval;
  if (!is_reload)
    EM_problem->start_time = start_time;
  else
    start_time = EM_problem->start_time;
  EM_problem->n_time_steps = (finish_time - start_time) / dt_integration;

  dtC = MIN(EM_problem->dz, EM_problem->dr) / c_mks / sqrt(2.0);
  if (dtC < EM_problem->dti)
    {
      if (auto_max_dt)
        {
          printf("warning: dt_integration (%e) exceeds max stable dt (%e); using max dt\n", dt_integration, dtC);
          EM_problem->dti = dtC;
          EM_problem->n_time_steps = (finish_time - start_time) / EM_problem->dti;
        }
      else
        {
          printf("error: the time step must be less than %es for stable integration.\n", dtC);
          exit(1);
        }
    }

  Jm1 = 1;
  Jm2 = 0;
  if (J_filter_multiplier)
    {
      Jm2 = 1 - J_filter_multiplier;
      Jm1 = J_filter_multiplier / Jm2;
    }

  printf("starting integration from time %e\n", EM_problem->start_time);
  if (terminate_on_total_loss)
    printf("integration will stop when all particles are lost.\n");
  printf("integration will consist of up to %ld steps spaced by %e\n\n",
         EM_problem->n_time_steps, EM_problem->dti);
  fflush(stdout);

  calculate_areas(EM_problem);

  beam_seen = 0;
  compute_imposed_Efield(EM_problem, 1);
  add_poisson_correction(EM_problem, beam, EM_problem->modes & FL_SPACE_CHARGE);
  for (icat = 0; icat < cathodes; icat++)
    set_current_density(cathode + icat, EM_problem);
  for (icat = 0; icat < nEmitters; icat++)
    setEmitterCurrentDensity(emitter + icat, EM_problem);
  for (; EM_problem->time_step < EM_problem->n_time_steps; EM_problem->time_step++)
    {
      if (EM_problem->imposedEFieldRampTime > 0)
        compute_imposed_Efield(EM_problem, 0);
      EM_problem->time = EM_problem->time_step * EM_problem->dti + EM_problem->start_time;
      if (terminate_on_total_loss)
        {
          if (beam_seen && !beam->np)
            {
              printf("\nsimulation terminating due to total loss of particles.\n");
              break;
            }
          if (beam->np && !beam_seen)
            beam_seen = 1;
        }
      sample_fields(sample_defn, n_sample_defns, EM_problem);
      save_fields(field_saving_defn, EM_problem, 0);
      output_fields(field_output_defn, EM_problem);
      take_beam_snapshot(snapshot_defn, beam, EM_problem);

      if (EM_problem->status_interval > 0 && (EM_problem->status_interval == 0 || EM_problem->time_step % EM_problem->status_interval == 0))
        do_status_printouts(EM_problem, beam, cathode,
                            check_divergence, EM_problem->modes & FL_SPACE_CHARGE,
                            status_output ? &SDDSstatus : NULL, &outputRow);

      if (nSecondaryDefns)
        emitSecondaryElectrons(beam, secondaryEmission, nSecondaryDefns, EM_problem);

      for (icat = 0; icat < cathodes; icat++)
        emit_electrons_from_cathode(beam, cathode + icat, EM_problem, icat);
      for (icat = 0; icat < nEmitters; icat++)
        emit_electrons_from_emitter(beam, emitter + icat, EM_problem);
      emit_test_particles(beam, test_beam_defn, EM_problem);

      if (!J_filter_multiplier)
        zero_current_arrays(EM_problem);
      else
        multiply_current_arrays(EM_problem, Jm1);
      if (EM_problem->modes & FL_SPACE_CHARGE)
        add_particle_currents(EM_problem, beam);
      add_antenna_excitation(EM_problem, antenna, n_antennae, EM_problem->time + EM_problem->dti / 2);
      if (resistor->n_resistors)
        add_resistor_currents(EM_problem, resistor);
      if (J_filter_multiplier)
        multiply_current_arrays(EM_problem, Jm2);
      if (smoothing_parameter)
        smooth_currents(EM_problem, smoothing_parameter);
      if (SG_smoothing_halfwidth)
        smooth_currents_SG(EM_problem, SG_smoothing_halfwidth);

      add_poisson_correction(EM_problem, beam, EM_problem->modes & FL_SPACE_CHARGE);
      advance_field_solutions(EM_problem);
      advance_particles(beam, EM_problem, test_beam_defn, screen_defn, n_screen_defns,
                        lost_particles ? &SDDSlost : NULL, z_forces_start,
                        secondaryEmission, nSecondaryDefns);
      translate_problem(beam, EM_problem, antenna, n_antennae, resistor);
      for (icat = 0; icat < cathodes; icat++)
        set_current_density(cathode + icat, EM_problem);
      for (icat = 0; icat < nEmitters; icat++)
        setEmitterCurrentDensity(emitter + icat, EM_problem);
    }
  EM_problem->time = EM_problem->time_step * EM_problem->dti + EM_problem->start_time;
  sample_fields(sample_defn, n_sample_defns, EM_problem);
  save_fields(field_saving_defn, EM_problem, 1);
  take_beam_snapshot(snapshot_defn, beam, EM_problem);
  if (EM_problem->status_interval > 0 && (EM_problem->status_interval == 0 || EM_problem->time_step % EM_problem->status_interval == 0))
    do_status_printouts(EM_problem, beam, cathode, check_divergence, EM_problem->modes & FL_SPACE_CHARGE,
                        status_output ? &SDDSstatus : NULL, &outputRow);

  if (lost_particles && !SDDS_Terminate(&SDDSlost))
    {
      SDDS_SetError("Problem terminate lost particle file");
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
    }
}

void zero_current_arrays(FIELDS *EM_problem)
{
  long iz, ir;

  if (EM_problem->modes & FL_TM_FIELDS)
    {
      if (!EM_problem->Jz)
        bomb("Jz is unexpectedly NULL in zero_current_arrays", NULL);
      if (!EM_problem->Jr)
        bomb("Jr is unexpectedly NULL in zero_current_arrays", NULL);
      for (iz = 0; iz <= EM_problem->nz; iz++)
        {
          for (ir = 0; ir <= EM_problem->nr; ir++)
            {
              EM_problem->Jz[iz][ir] = EM_problem->Jr[iz][ir] = 0;
            }
        }
    }

  if (EM_problem->modes & FL_TE_FIELDS)
    {
      if (!EM_problem->Jphi)
        bomb("Jphi is unexpectedly NULL in zero_current_arrays", NULL);
      for (iz = 0; iz <= EM_problem->nz; iz++)
        {
          for (ir = 0; ir <= EM_problem->nr; ir++)
            {
              EM_problem->Jphi[iz][ir] = 0;
            }
        }
    }
}

void multiply_current_arrays(FIELDS *EM_problem, double factor)
{
  long iz, ir;

  if (EM_problem->modes & FL_TM_FIELDS)
    {
      if (!EM_problem->Jz)
        bomb("Jz is unexpectedly NULL in multiply_current_arrays", NULL);
      if (!EM_problem->Jr)
        bomb("Jr is unexpectedly NULL in multiply_current_arrays", NULL);
      for (iz = 0; iz <= EM_problem->nz; iz++)
        {
          for (ir = 0; ir <= EM_problem->nr; ir++)
            {
              EM_problem->Jz[iz][ir] *= factor;
              EM_problem->Jr[iz][ir] *= factor;
            }
        }
    }

  if (EM_problem->modes & FL_TE_FIELDS)
    {
      for (iz = 0; iz <= EM_problem->nz; iz++)
        {
          for (ir = 0; ir <= EM_problem->nr; ir++)
            {
              EM_problem->Jphi[iz][ir] *= factor;
            }
        }
      if (!EM_problem->Jphi)
        bomb("Jphi is unexpectedly NULL in multiply_current_arrays", NULL);
    }
}

void advance_field_solutions(FIELDS *EM_problem)
{
  long iz, ir, nz, nr;
  double dr, dz, dti, dB;
  double **Ez, **Er, **Bphi, **Bphi2, **Jz, **Jr;
  double **Ephi, **Jphi, **Br, **Br2, **Bz, **Bz2;
  short **metal_flags;
  double rp, rm, factor, mu0, eps0, cp, cm;

  nz = EM_problem->nz;
  nr = EM_problem->nr;
  dr = EM_problem->dr;
  dz = EM_problem->dz;
  dti = EM_problem->dti;
  eps0 = epsilon_o;
  mu0 = mu_o;

  Ez = EM_problem->Ez;
  Jz = EM_problem->Jz;
  Er = EM_problem->Er;
  Jr = EM_problem->Jr;
  Bphi = EM_problem->Bphi;
  Bphi2 = EM_problem->Bphi2;
  Ephi = EM_problem->Ephi;
  Jphi = EM_problem->Jphi;
  Bz = EM_problem->Bz;
  Bz2 = EM_problem->Bz2;
  Br = EM_problem->Br;
  Br2 = EM_problem->Br2;
  metal_flags = EM_problem->metal_flags;

  if (EM_problem->modes & FL_TM_FIELDS)
    {
      /* advance TM fields */
      if (!EM_problem->Ez)
        bomb("Ez is unexpectedly NULL in advance_field_solutions", NULL);
      if (!EM_problem->Jz)
        bomb("Jz is unexpectedly NULL in advance_field_solutions", NULL);
      if (!EM_problem->Er)
        bomb("Er is unexpectedly NULL in advance_field_solutions", NULL);
      if (!EM_problem->Jr)
        bomb("Jr is unexpectedly NULL in advance_field_solutions", NULL);
      if (!EM_problem->Bphi)
        bomb("Bphi is unexpectedly NULL in advance_field_solutions", NULL);
      if (!EM_problem->Bphi2)
        bomb("Bphi2 is unexpectedly NULL in advance_field_solutions", NULL);
      /* -- Ez */
      factor = 1 / (eps0 * mu0);
      for (ir = 0; ir < nr; ir++)
        {
          rp = (ir + 0.5) * dr;
          rm = (ir - 0.5) * dr;
          if (ir != 0)
            {
              /* should store these for reuse */
              cp = 2 * factor / (sqr(rp) - sqr(rm));
              cm = rm * cp;
              cp = rp * cp;
            }
          else
            {
              cp = 4 * factor / dr;
              cm = 0;
            }
          for (iz = 1; iz < nz; iz++)
            {
              if (!(metal_flags[iz][ir] & FL_EZ_ZERO))
                Ez[iz][ir] += dti * (-Jz[iz][ir] / eps0 + (cp * Bphi[iz][ir + 1] - cm * Bphi[iz][ir]));
              if (metal_flags[iz][ir] & FL_LAST_ACTIVE_PT)
                break;
            }
          if (EM_problem->left_bc == NEUMANN)
            Ez[0][ir] = Ez[1][ir];
          else
            Ez[0][ir] = -Ez[1][ir];
          if (EM_problem->right_bc == NEUMANN)
            Ez[nz][ir] = Ez[nz - 1][ir];
          else
            Ez[nz][ir] = -Ez[nz - 1][ir];
        }
      /* Neumann B.C. on upper boundary */
      for (iz = 0; iz <= nz; iz++)
        Ez[iz][nr - 1] = 0;

      /* -- Er */
      factor = 1 / (mu0 * eps0 * dz);
      for (ir = 1; ir < nr; ir++)
        {
          for (iz = 0; iz < nz; iz++)
            {
              if (!(metal_flags[iz][ir] & FL_ER_ZERO))
                Er[iz][ir] += dti * (-Jr[iz][ir] / eps0 - factor * (Bphi[iz + 1][ir] - Bphi[iz][ir]));
              if (metal_flags[iz][ir] & FL_LAST_ACTIVE_PT)
                break;
            }
        }
      for (iz = 0; iz < nz; iz++)
        Er[iz][0] = -Er[iz][1]; /* Dirichlet conditions on lower boundary */
      if (EM_problem->left_bc == NEUMANN)
        for (ir = 0; ir <= nr; ir++)
          Er[0][ir] = 0;
      if (EM_problem->right_bc == NEUMANN)
        for (ir = 0; ir <= nr; ir++)
          Er[nz - 1][ir] = 0;

      /* advance magnetic field values */
      /* -- Bphi */
      factor = (EM_problem->rmax - EM_problem->dr / 2) / (EM_problem->rmax + EM_problem->dr / 2);
      for (ir = 1; ir < nr; ir++)
        {
          for (iz = 1; iz < nz; iz++)
            {
              if (!(metal_flags[iz][ir] & FL_BPHI_ZERO))
                {
                  dB = dti * ((Ez[iz][ir] - Ez[iz][ir - 1]) / dr -
                              (Er[iz][ir] - Er[iz - 1][ir]) / dz);
                  Bphi2[iz][ir] = Bphi[iz][ir] + dB / 2;
                  Bphi[iz][ir] += dB;
                }
              if (metal_flags[iz][ir] & FL_LAST_ACTIVE_PT)
                break;
            }
        }
      for (iz = 1; iz < nz; iz++)
        {
          /* Dirichlet boundary at bottom */
          Bphi[iz][0] = -Bphi[iz][1];
          Bphi2[iz][0] = -Bphi2[iz][1];
          /* Neumann boundary at top */
          Bphi[iz][nr - 1] = Bphi[iz][nr - 2] * factor;
          Bphi2[iz][nr - 1] = Bphi2[iz][nr - 2] * factor;
        }
      if (EM_problem->left_bc == NEUMANN)
        for (ir = 0; ir <= nr; ir++)
          {
            Bphi[0][ir] = Bphi[1][ir];
            Bphi2[0][ir] = Bphi2[1][ir];
          }
      else
        for (ir = 0; ir <= nr; ir++)
          {
            Bphi[0][ir] = -Bphi[1][ir];
            Bphi2[0][ir] = -Bphi2[1][ir];
          }
      if (EM_problem->right_bc == NEUMANN)
        for (ir = 0; ir <= nr; ir++)
          {
            Bphi[nz][ir] = Bphi[nz - 1][ir];
            Bphi2[nz][ir] = Bphi2[nz - 1][ir];
          }
      else
        for (ir = 0; ir <= nr; ir++)
          {
            Bphi[nz][ir] = -Bphi[nz - 1][ir];
            Bphi2[nz][ir] = -Bphi2[nz - 1][ir];
          }
    }
  if (EM_problem->modes & FL_TE_FIELDS)
    {
      /* advance TE fields */
      if (!EM_problem->Ephi)
        bomb("Ephi is unexpectedly NULL in advance_field_solutions", NULL);
      if (!EM_problem->Jphi)
        bomb("Jphi is unexpectedly NULL in advance_field_solutions", NULL);
      if (!EM_problem->Br)
        bomb("Br is unexpectedly NULL in advance_field_solutions", NULL);
      if (!EM_problem->Br2)
        bomb("Br2 is unexpectedly NULL in advance_field_solutions", NULL);
      if (!EM_problem->Bz)
        bomb("Bz is unexpectedly NULL in advance_field_solutions", NULL);
      if (!EM_problem->Bz2)
        bomb("Bz2 is unexpectedly NULL in advance_field_solutions", NULL);
      /* -- Ephi */
      factor = 1. / (eps0 * mu0);
      for (ir = 1; ir < nr; ir++)
        for (iz = 0; iz < nz; iz++)
          {
            if (!(metal_flags[iz][ir] & FL_EPHI_ZERO))
              Ephi[iz][ir] += dti * (-Jphi[iz][ir] / eps0 +
                                     factor * ((Br[iz + 1][ir] - Br[iz][ir]) / dz - (Bz[iz][ir + 1] - Bz[iz][ir]) / dr));
            if (metal_flags[iz][ir] & FL_LAST_ACTIVE_PT)
              break;
          }
      /* Dirichlet B.C. along lower boundary */
      for (iz = 0; iz <= nz; iz++)
        Ephi[iz][0] = 0;
      if (EM_problem->right_bc == NEUMANN)
        for (ir = 0; ir <= nr; ir++)
          Ephi[nz - 1][ir] = 0;
      if (EM_problem->left_bc == NEUMANN)
        for (ir = 0; ir <= nr; ir++)
          Ephi[0][ir] = 0;
      if (EM_problem->upper_bc == NEUMANN)
        for (iz = 0; iz <= nz; iz++)
          Ephi[iz][nr - 1] = 0;

      /* -- Br */
      for (ir = 0; ir < nr; ir++)
        for (iz = 1; iz < nz; iz++)
          {
            if (!(metal_flags[iz][ir] & FL_BR_ZERO))
              {
                dB = dti * (Ephi[iz][ir] - Ephi[iz - 1][ir]) / dz;
                Br2[iz][ir] = Br[iz][ir] + dB / 2;
                Br[iz][ir] += dB;
              }
            if (metal_flags[iz][ir] & FL_LAST_ACTIVE_PT)
              break;
          }
      if (EM_problem->right_bc == NEUMANN)
        for (ir = 0; ir <= nr; ir++)
          {
            Br[nz][ir] = Br[nz - 1][ir];
            Br2[nz][ir] = Br2[nz - 1][ir];
          }
      else
        for (ir = 0; ir <= nr; ir++)
          {
            Br[nz][ir] = -Br[nz - 1][ir];
            Br2[nz][ir] = -Br2[nz - 1][ir];
          }
      if (EM_problem->left_bc == NEUMANN)
        for (ir = 0; ir <= nr; ir++)
          {
            Br[0][ir] = Br[1][ir];
            Br2[0][ir] = Br2[1][ir];
          }
      else
        for (ir = 0; ir <= nr; ir++)
          {
            Br[0][ir] = -Br[1][ir];
            Br2[0][ir] = -Br2[1][ir];
          }

      /* -- Bz */
      for (ir = 1; ir < nr; ir++)
        {
          /* should store these factors for later use */
          factor = -2 * dti / (dr * (isqr(ir + 1) - isqr(ir)));
          for (iz = 0; iz < nz; iz++)
            {
              if (!(metal_flags[iz][ir] & FL_BZ_ZERO))
                {
                  dB = factor * ((ir + 1) * Ephi[iz][ir] - ir * Ephi[iz][ir - 1]);
                  Bz2[iz][ir] = Bz[iz][ir] + dB / 2;
                  Bz[iz][ir] += dB;
                }
              if (metal_flags[iz][ir] & FL_LAST_ACTIVE_PT)
                break;
            }
        }
      /* Dirichlet B.C. along lower boundary */
      for (iz = 0; iz <= nz; iz++)
        {
          Bz[iz][0] = Bz[iz][1];
          Bz2[iz][0] = Bz2[iz][1];
        }
      if (EM_problem->right_bc == NEUMANN)
        for (ir = 0; ir <= nr; ir++)
          Bz[nz - 1][ir] = Bz2[nz - 1][ir] = 0;
      if (EM_problem->left_bc == NEUMANN)
        for (ir = 0; ir <= nr; ir++)
          Bz[0][ir] = Bz2[0][ir] = 0;
    }
}

void smooth_currents(
                     FIELDS *EM_problem,
                     double smoothing_parameter)
{
  long iz, ir, nz, nr;
  double **Jz, **Jr, **Jphi;
  static double *buffer = NULL;
  /*short **metal_flags;*/
  double c1, c2;

  c1 = (1 - smoothing_parameter);
  c2 = smoothing_parameter / 2;

  nz = EM_problem->nz;
  nr = EM_problem->nr;

  Jz = EM_problem->Jz;
  Jr = EM_problem->Jr;
  Jphi = EM_problem->Jphi;
  /*metal_flags = EM_problem->metal_flags;*/

  if (!buffer)
    buffer = (double *)tmalloc(sizeof(*buffer) * MAX(nz + 1, nr + 1));

  if (EM_problem->modes & FL_TM_FIELDS)
    {
      if (!EM_problem->Jz)
        bomb("Jz array is unexpectedly NULL in smooth_currents", NULL);
      if (!EM_problem->Jr)
        bomb("Jr array is unexpectedly NULL in smooth_currents", NULL);
      /* smooth Jz */
      for (ir = 0; ir < nr; ir++)
        {
          /* smooth z variation for constant r */
          buffer[1] = (c1 + c2) * Jz[1][ir] + c2 * Jz[2][ir];
          for (iz = 2; iz < nz - 1; iz++)
            buffer[iz] = c1 * Jz[iz][ir] + c2 * (Jz[iz - 1][ir] + Jz[iz + 1][ir]);
          buffer[nz - 1] = (c1 + c2) * Jz[nz - 1][ir] + c2 * Jz[nz - 2][ir];
          for (iz = 1; iz < nz; iz++)
            Jz[iz][ir] = buffer[iz];
        }

      /* smooth Jr */
      for (ir = 1; ir < nr; ir++)
        {
          /* smooth z variation for constant r */
          buffer[0] = (c1 + c2) * Jr[0][ir] + c2 * Jr[1][ir];
          for (iz = 1; iz < nz - 1; iz++)
            buffer[iz] = c1 * Jr[iz][ir] + c2 * (Jr[iz - 1][ir] + Jr[iz + 1][ir]);
          buffer[nz - 1] = (c1 + c2) * Jr[nz - 1][ir] + c2 * Jr[nz - 2][ir];
          for (iz = 0; iz < nz; iz++)
            Jr[iz][ir] = buffer[iz];
        }
    }

  if (EM_problem->modes & FL_TE_FIELDS)
    {
      if (!EM_problem->Jphi)
        bomb("Jphi array is unexpectedly NULL in smooth_currents", NULL);
      /* smooth Jphi */
      for (ir = 0; ir < nr; ir++)
        {
          /* smooth z variation for constant r */
          buffer[0] = (c1 + c2) * Jphi[0][ir] + c2 * Jphi[1][ir];
          for (iz = 1; iz < nz - 1; iz++)
            buffer[iz] = c1 * Jphi[iz][ir] + c2 * (Jphi[iz - 1][ir] + Jphi[iz + 1][ir]);
          buffer[nz - 1] = (c1 + c2) * Jphi[nz - 1][ir] + c2 * Jphi[nz - 2][ir];
          for (iz = 0; iz < nz; iz++)
            Jphi[iz][ir] = buffer[iz];
        }
    }
}

void smooth_currents_SG(
                        FIELDS *EM_problem,
                        long SG_HW)
{
  long iz, ir, nz, nr;
  double **Jz, **Jr, **Jphi, sum;
  static double *buffer = NULL;
  /*short **metal_flags;*/

  nz = EM_problem->nz;
  nr = EM_problem->nr;

  Jz = EM_problem->Jz;
  Jr = EM_problem->Jr;
  Jphi = EM_problem->Jphi;
  /*metal_flags = EM_problem->metal_flags;*/

  if (!buffer)
    buffer = (double *)tmalloc(sizeof(*buffer) * MAX(nz + 1, nr + 1));

  if (EM_problem->modes & FL_TM_FIELDS)
    {
      if (!EM_problem->Jz)
        bomb("Jz array is unexpectedly NULL in smooth_currents", NULL);
      if (!EM_problem->Jr)
        bomb("Jr array is unexpectedly NULL in smooth_currents", NULL);
      /* smooth Jz */
      for (ir = 0; ir < nr; ir++)
        {
          sum = 0;
          for (iz = 1; iz < nz; iz++)
            sum += fabs(buffer[iz - 1] = Jz[iz][ir]);
          if (sum)
            {
              SavitzkyGolaySmooth(buffer, nz - 1, 1, SG_HW, SG_HW, 0);
              for (iz = 1; iz < nz; iz++)
                Jz[iz][ir] = buffer[iz - 1];
            }
        }

      /* smooth Jr */
      for (ir = 1; ir < nr; ir++)
        {
          for (iz = sum = 0; iz < nz; iz++)
            sum += fabs(buffer[iz] = Jr[iz][ir]);
          if (sum)
            {
              SavitzkyGolaySmooth(buffer, nz, 1, SG_HW, SG_HW, 0);
              for (iz = 0; iz < nz; iz++)
                Jr[iz][ir] = buffer[iz];
            }
        }
    }

  if (EM_problem->modes & FL_TE_FIELDS)
    {
      if (!EM_problem->Jphi)
        bomb("Jphi array is unexpectedly NULL in smooth_currents", NULL);
      /* smooth Jphi */
      for (ir = 0; ir < nr; ir++)
        {
          for (iz = sum = 0; iz < nz; iz++)
            sum += fabs(buffer[iz] = Jphi[iz][ir]);
          if (sum)
            {
              SavitzkyGolaySmooth(buffer, nz, 1, SG_HW, SG_HW, 0);
              for (iz = 0; iz < nz; iz++)
                Jphi[iz][ir] = buffer[iz];
            }
        }
    }
}

void calculate_areas(FIELDS *EM_problem)
{
  long ir;
  double dr2;

  dr2 = sqr(EM_problem->dr);
  /* calculate areas for annuli */
  for (ir = 0; ir < EM_problem->nr + 1; ir++)
    {
      if (ir)
        EM_problem->area_c[ir] = PIx2 * ir * dr2;
      else
        EM_problem->area_c[ir] = PI * dr2 / 4;
    }
  EM_problem->area_o[0] = 0;
  for (ir = 1; ir < EM_problem->nr + 1; ir++)
    EM_problem->area_o[ir] = PI * (2 * ir - 1) * dr2;
}
