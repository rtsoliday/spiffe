/* file: semaphores.nl
 * purpose: namelist for creation of semaphore files indicating run status.
 *          Modeled on the semaphores namelist in elegant.
 */
#include "namelist.h"

#namelist semaphores
    STRING started = "%s.started";
    STRING done = "%s.done";
    STRING failed = "%s.failed";
#end
