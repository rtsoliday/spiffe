/* file: load_fields.nl
 * purpose: command for loading saved fields
 *
 * Michael Borland, 1992 
 */
#include "namelist.h"

#namelist load_fields
    STRING filename = NULL;
    double Ez_peak;
    double factor;
    double time_threshold = 0;
    long overlay = 0;
#end


