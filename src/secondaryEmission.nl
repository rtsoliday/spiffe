/* file: secondary_emission.nl
 * purpose: namelist for setting up secondary emission.
 *
 * Michael Borland, 2005
 */
#include "namelist.h"

#namelist secondary_emission
        STRING input_file = NULL;
        STRING kinetic_energy_column = NULL;
        STRING yield_column = NULL;
        double yield_limit = 0;
        double emitted_momentum = 0;
        long verbosity = 1;
        long material_id = 1;
        STRING log_file = NULL;
#end
