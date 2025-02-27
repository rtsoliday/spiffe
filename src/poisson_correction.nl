/* file: poisson_correction.nl
 * purpose: setup for poisson correction
 *
 * Michael Borland, 1992
 */
#include "namelist.h"

#namelist poisson_correction
    double start_time;
    long step_interval;
    double accuracy;
    double error_charge_threshold;
    long maximum_iterations = 0;
    long verbosity;
    double test_charge;
    double z_test_charge;
    double r_test_charge;
    STRING guess_type;
#end

