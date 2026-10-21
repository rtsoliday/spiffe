/* file: rpnCommands.nl
 * purpose: namelists for rpn expression evaluation, rpn data loading,
 *          and macro output.  Modeled on elegant's run_rpnexpr.nl and
 *          the macro_output namelist in elegant.nl.
 */
#include "namelist.h"

#namelist rpn_expression
    STRING expression = NULL;
#end

#namelist rpn_load
    STRING filename = NULL;
    STRING match_column = NULL;
    STRING match_column_value = NULL;
    long matching_row_number = -1;
    STRING match_parameter = NULL;
    STRING match_parameter_value = NULL;
    long use_row = -1;
    long use_page = -1;
    long load_parameters = 0;
    STRING tag = NULL;
#end

#namelist macro_output,struct
    STRING filename = NULL;
    STRING mode = "parameters";
#end
