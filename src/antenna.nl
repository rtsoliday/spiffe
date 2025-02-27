/* file: antenna.nl
 * purpose: namelist for antenna setup
 *
 * Michael Borland, 1992
 */
#include "namelist.h"

#namelist define_antenna
    double start;
    double end;
    double position;
    STRING direction;     /* "z" ("r") implies r0=r1 (z0=z1) */
    double current;       /* in Amperes */
    double frequency;     /* in Hertz */
    double phase;         /* in radians */
    STRING waveform;      /* mpl-format file giving the excitation waveform, W(t) */
    double time_offset;   /* antenna is driven with current*W(t-time_offset)*sin(2*pi*frequency*t+phase) */
#end
