/* file: field_output.nl
 * purpose: setup for output of fields
 *
 * Michael Borland, 1992 
 */
#include "namelist.h"

#namelist define_field_output
    STRING filename;
    double time_interval;
    double start_time;
    long z_interval;
    long r_interval;
    long exclude_imposed_field;
    long separate_imposed_field;
#end



