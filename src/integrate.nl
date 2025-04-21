/* file: integration.nl
 * purpose: setup for integration
 *
 * Michael Borland, 1992
 */
#include "namelist.h"

#namelist integrate
    double dt_integration;   /* dt for integration */
    double start_time;       /* time at which integration starts */
    double finish_time;      /* number of steps is (finish_time-start_time)/dt_integration */
    long status_interval;    /* number of steps between status reports */
    long space_charge;       /* logical--if non-zero, currents due to particles are included */
    long check_divergence;   /* logical--if non-zero, div.E is computed and reported in the status reports */
    double smoothing_parameter;  /* A[i] -> (1-s.p.)*A[i] + s.p./2*(A[i-1]+A[i+1]) */
    double J_filter_multiplier;  
    long SG_smoothing_halfwidth;  /* another way to smooth J */
    long terminate_on_total_loss;  /* if non-zero, will terminate when all particles are lost */
    double z_forces_start = -1;
    STRING status_output = NULL;
    STRING lost_particles = NULL;
    double imposed_Efield_ramp_time = 0;
    double imposed_Efield_flat_top_time = 1;
    long auto_max_dt = 0;    /* if non-zero, use maximum stable dt when dt_integration too large */
#end
