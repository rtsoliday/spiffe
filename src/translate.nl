/* file: translate.nl
 * purpose: setup for z translation of mesh and particles
 *
 * Michael Borland, 1992
 */
#include "namelist.h"

#namelist translate
    double z_trigger;        /* translation begins when a particle passes this z */
    double z_lower;          /* lower limit of region to include in simulation after translation */
    double z_upper;          /* upper limit of same region */
#end
