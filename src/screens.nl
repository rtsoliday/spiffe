/* file: screens.nl
 * purpose: setup for beam-imaging screens
 *
 * Michael Borland, 1992 
 */
#include "namelist.h"

#namelist define_screen
    STRING filename;
    STRING template;
    double z_position;
    double delta_z;
    long number_of_screens;
    double start_time;
    STRING direction;
#end


