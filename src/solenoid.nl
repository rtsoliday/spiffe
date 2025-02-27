/* file: solenoid.nl
 * purpose: define a solenoid coil
 *
 * Michael Borland, 1995
 */
#include "namelist.h"

#namelist define_solenoid
        double radius;
        double evaluation_radius_limit;
        double z_start;
        double z_end;
        double current;
        double Bz_peak;
        long turns;
        long symmetry;
        STRING field_output = NULL;
        long bucking = 0;
        double z_buck = 0;
#end

