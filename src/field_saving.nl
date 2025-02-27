/* file: field_saving.nl
 * purpose: setup for saving of fields
 *
 * Michael Borland, 1992 
 */
#include "namelist.h"

#namelist define_field_saving
    STRING filename = NULL;
    double time_interval = -1;
    double start_time = 0;
    long save_before_exiting = 0;
    double z_Ez_trigger = -1;
#end


