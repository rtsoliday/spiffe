/* file: field_sampling.nl
 * purpose: setup for sampling of fields
 *
 * Michael Borland, 1992 
 */
#include "namelist.h"

#namelist define_field_sampling
    STRING filename;
    STRING component;
    STRING direction;
    double min_coord;
    double max_coord;
    double position;
    double time_interval;
    double start_time;
    long time_sequence;
#end


