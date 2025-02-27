/* file: emitter.nl
 * purpose: namelist for defining emitter region
 *
 * Michael Borland, 2017
 */
#include "namelist.h"

#namelist define_emitter
    long material_id = -1;
    double current_density;
    double temperature;
    long determine_temperature;
    double work_function ; /* in eV */
    long add_thermal_velocities = 0;
    long field_emission = 0;
    double field_emission_beta = 1;
    double electrons_per_macroparticle;
    double start_time;
    double stop_time;
    double number_per_step;
    double initial_p;
    double offset_factor = 1e-2;
    int random_number_seed;
    int32_t halton_radix_dt = 0;
    int32_t halton_radix_r = 0;
    int32_t halton_radix_z = 0;
    STRING profile = NULL;
    STRING profile_time_name = NULL;
    STRING profile_factor_name = NULL;
    STRING emission_log = NULL;
#end
