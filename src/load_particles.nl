/* file: load_particles.nl
 * purpose: command for loading particles from an external file
 *
 * Michael Borland, 1992 
 */
#include "namelist.h"

#namelist load_particles
        STRING filename = NULL;
        long sample_interval = 1;
        double stiffness = 1;
#end


