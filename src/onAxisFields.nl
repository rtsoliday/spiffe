/* file: onAxisField.nl
 * purpose: command for defining on-axis fields for pushing particles
 *
 * Michael Borland, 2003
 */
#include "namelist.h"

#namelist add_on_axis_fields
    STRING filename = NULL;
    STRING z_name = NULL;
    STRING Ez_name = NULL;
    double Ez_peak = 0;
    double frequency = 0;
    double phase = 0;
    double z_offset = 0;
    long expansion_order = 3;
    STRING fields_used = NULL;
#end


