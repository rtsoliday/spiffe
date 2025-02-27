/* file: test_beam.nl
 * purpose: namelist for defining test beam
 *
 * Michael Borland, 1992
 */
#include "namelist.h"

#namelist define_test_beam
    double start_time;
    double stop_time;
    double initial_z;
    double initial_r;
    double initial_pz;
    double initial_pr;
    double initial_rpphi;
    double z_force;
    double r_force;
    double z_start_deceleration;
    double r_start_deceleration;
    double electrons_per_macroparticle;
#end
