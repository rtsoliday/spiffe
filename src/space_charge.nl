/* file: space_charge.nl
 * purpose: setup for space charge--this is an obsolete namelist (flags are now in geometry.nl)
 *
 * Michael Borland, 1992
 */
#include "namelist.h"

#namelist space_charge
    long radial_sharing = 0;
    long longitudinal_sharing = 1;
#end

