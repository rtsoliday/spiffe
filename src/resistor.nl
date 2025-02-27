/* file: resistor.nl
 * purpose: namelist for resistor setup
 *
 * Michael Borland, 1992
 */
#include "namelist.h"

#namelist define_resistor
    double start;
    double end;
    double position;
    STRING direction;     /* "z" ("r") implies r0=r1 (z0=z1) */
    double conductivity;  /* in 1/(m*ohm) */
#end
