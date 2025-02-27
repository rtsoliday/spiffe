/* file: cathode.nl
 * purpose: namelist for defining cathode
 *
 * Michael Borland, 1992
 */
#include "namelist.h"

#namelist define_cathode
    double z_position;
    double inner_radius;
    double outer_radius;
    double current_density;
    double temperature;
    double work_function ; /* in eV */
    long field_emission = 0;
    long add_thermal_velocities = 0;
    double field_emission_beta = 1;
    long determine_temperature;
    double electrons_per_macroparticle;
    double start_time;
    double stop_time;
    long autophase;
    double time_offset;
    double number_per_step;
    double initial_pz;
    double initial_omega;
    double focal_length;
    double stiffness;
    long discretize_radii;
    int random_number_seed;
    long distribution_correction_interval;
    long spread_over_dt;
    long zoned_emission;
    int32_t halton_radix_dt = 0;
    int32_t halton_radix_r = 0;
    STRING profile = NULL;
    STRING profile_time_name = NULL;
    STRING profile_factor_name = NULL;
    STRING emission_log = NULL;
#end
