/*************************************************************************\
 * Copyright (c) 2002 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2002 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution.
\*************************************************************************/

/* file: rpnCommands.c
 * purpose: rpn_expression, rpn_load, and macro_output commands.
 *          Ported from elegant's run_rpnexpr.c and elegant.c.
 */
#include "spiffe.h"
#include "match_string.h"
#include "rpn.h"
#include "rpnCommands.h"
#include <ctype.h>

void run_rpn_expression(NAMELIST_TEXT *nl_text)
{
  expression = NULL;

  process_namelist(&rpn_expression, nl_text);
  print_namelist(stdout, &rpn_expression);

  if (expression)
    {
      rpn(expression);
      if (rpn_check_error())
        exit(1);
    }
}

void run_rpn_load(NAMELIST_TEXT *nl_text, FIELDS *EM_problem)
{
  SDDS_DATASET SDDSin;
  long code, foundPage, iColumn, matchRow, rows, iParameter;
  int32_t columns, parameters;
  char *parameterValue = NULL;
  double *data, data1;
  char **columnName, **matchColumnData, *memName = NULL;
  char **parameterName;

  filename = match_column = match_column_value = NULL;
  match_parameter = match_parameter_value = tag = NULL;
  matching_row_number = use_row = use_page = -1;
  load_parameters = 0;

  process_namelist(&rpn_load, nl_text);
  print_namelist(stdout, &rpn_load);

  if (!tag || strlen(tag) == 0)
    fprintf(stderr, "Warning: rpn_load tag is blank, which may cause confusion with internal variables.\n");

  if (match_column && strlen(match_column))
    {
      if (use_row != -1)
        bomb("rpn_load: you asked to match a column and also gave use_row", NULL);
      if (!match_column_value || !strlen(match_column_value))
        bomb("rpn_load: you must give match_column_value with match_column", NULL);
    }
  if (match_parameter && strlen(match_parameter))
    {
      if (use_page != -1)
        bomb("rpn_load: you asked to match a parameter and also gave use_page", NULL);
      if (!match_parameter_value || !strlen(match_parameter_value))
        bomb("rpn_load: you must give match_parameter_value with match_parameter", NULL);
    }

  if (!filename || !strlen(filename))
    bomb("rpn_load: no filename given", NULL);

  filename = compose_filename(filename, EM_problem->rootname);

  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, filename))
    {
      fprintf(stderr, "Error: couldn't initialize SDDS input for %s\n", filename);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }

  foundPage = 0;
  while ((code = SDDS_ReadPage(&SDDSin)) > 0)
    {
      if (use_page > 0)
        {
          if (code == use_page)
            {
              foundPage = 1;
              break;
            }
          continue;
        }
      if (match_parameter && strlen(match_parameter))
        {
          if (!(parameterValue = SDDS_GetParameterAsString(&SDDSin, match_parameter, NULL)))
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors | SDDS_EXIT_PrintErrors);
          if (!wild_match(parameterValue, match_parameter_value))
            continue;
          foundPage = 1;
          break;
        }
      if (use_page == -1 && SDDS_CheckEndOfFile(&SDDSin) == 1)
        {
          foundPage = 1;
          break;
        }
    }

  if (!foundPage)
    bomb("rpn_load: no appropriate page found", NULL);

  if (!load_parameters)
    {
      if ((columnName = SDDS_GetColumnNames(&SDDSin, &columns)) == NULL)
        {
          fprintf(stderr, "Warning: rpn_load found no columns in %s\n", filename);
          return;
        }

      rows = SDDS_RowCount(&SDDSin);
      matchRow = rows - 1;
      if (use_row != -1)
        {
          if (use_row >= rows)
            {
              fprintf(stderr, "Error: number of rows in file (%ld) less than needed for use_row=%ld\n",
                      rows, use_row);
              exit(1);
            }
          matchRow = use_row;
        }

      if (match_column)
        {
          if (SDDS_GetNamedColumnType(&SDDSin, match_column) != SDDS_STRING)
            {
              fprintf(stderr, "Error: column %s nonexistent or not string type.\n", match_column);
              exit(1);
            }
          if (!(matchColumnData = SDDS_GetColumn(&SDDSin, match_column)))
            {
              fprintf(stderr, "Error: unable to get data for column %s\n", match_column);
              exit(1);
            }
          if (matching_row_number < 0)
            {
              /* use last match */
              for (matchRow = rows - 1; matchRow >= 0; matchRow--)
                if (wild_match(matchColumnData[matchRow], match_column_value))
                  break;
            }
          else
            {
              /* use nth match */
              for (matchRow = 0; matchRow < rows; matchRow++)
                if (wild_match(matchColumnData[matchRow], match_column_value) &&
                    matching_row_number-- == 0)
                  break;
            }

          if (matchRow < 0 || matchRow >= rows)
            {
              fprintf(stderr, "Error: unable to find match for %s in column %s\n",
                      match_column_value, match_column);
              exit(1);
            }
          SDDS_FreeStringArray(matchColumnData, rows);
        }

      for (iColumn = 0; iColumn < columns; iColumn++)
        {
          switch (SDDS_GetNamedColumnType(&SDDSin, columnName[iColumn]))
            {
            case SDDS_CHARACTER:
            case SDDS_STRING:
              break;
            default:
              if (!(data = SDDS_GetColumnInDoubles(&SDDSin, columnName[iColumn])))
                {
                  fprintf(stderr, "Error: unable to get data for column %s as numerical data.\n",
                          columnName[iColumn]);
                  exit(1);
                }
              if (!(memName = SDDS_Realloc(memName, sizeof(*memName) * ((tag ? strlen(tag) : 0) + strlen(columnName[iColumn]) + 2))))
                bomb("memory allocation failure creating memory name for loaded data", NULL);
              if (tag && strlen(tag))
                sprintf(memName, "%s.%s", tag, columnName[iColumn]);
              else
                sprintf(memName, "%s", columnName[iColumn]);
              rpn_store(data[matchRow], NULL, rpn_create_mem(memName, 0));
              printf("%le --> %s\n", data[matchRow], memName);
              free(columnName[iColumn]);
              free(data);
            }
        }
      if (memName)
        free(memName);
      if (columnName)
        free(columnName);
    }
  else
    {
      /* load data from parameters */
      if ((parameterName = SDDS_GetParameterNames(&SDDSin, &parameters)) == NULL)
        {
          fprintf(stderr, "Warning: rpn_load found no parameters in %s\n", filename);
          return;
        }

      for (iParameter = 0; iParameter < parameters; iParameter++)
        {
          switch (SDDS_GetNamedParameterType(&SDDSin, parameterName[iParameter]))
            {
            case SDDS_CHARACTER:
            case SDDS_STRING:
              break;
            default:
              if (!SDDS_GetParameterAsDouble(&SDDSin, parameterName[iParameter], &data1))
                {
                  fprintf(stderr, "Error: unable to get data for parameter %s as numerical data.\n",
                          parameterName[iParameter]);
                  exit(1);
                }
              if (!(memName = SDDS_Realloc(memName, sizeof(*memName) * ((tag ? strlen(tag) : 0) + strlen(parameterName[iParameter]) + 2))))
                bomb("memory allocation failure creating memory name for loaded data", NULL);
              if (tag && strlen(tag))
                sprintf(memName, "%s.%s", tag, parameterName[iParameter]);
              else
                sprintf(memName, "%s", parameterName[iParameter]);
              rpn_store(data1, NULL, rpn_create_mem(memName, 0));
              printf("%le --> %s\n", data1, memName);
              free(parameterName[iParameter]);
            }
        }
      if (memName)
        free(memName);
      if (parameterName)
        free(parameterName);
    }

  if (!SDDS_Terminate(&SDDSin))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
}

/* Determine the narrowest SDDS type that can represent a macro value string.
 * Copied from elegant's SDDS_IdentifyTypeOfData() so that macro_output files
 * from spiffe and elegant are typed identically.
 */
static long identifyTypeOfData(char *data0)
{
#define DATA_IS_INT 0x01UL
#define DATA_IS_FLOAT 0x02UL
  unsigned long flags = DATA_IS_INT | DATA_IS_FLOAT;
  short exponentialSeen = 0, periodSeen = 0;
  char *data;
  long longValue;
  double doubleValue;

  data = data0;
  if (*data == '+' || *data == '-')
    data++;

  while (*data && flags && !exponentialSeen)
    {
      if (*data <= '/' && !(*data == '.' && !periodSeen))
        {
          flags = 0;
          break;
        }
      if (*data >= ':' || *data == '.')
        {
          flags &= ~DATA_IS_INT;
          if (*data == '.' && flags & DATA_IS_FLOAT)
            {
              if (periodSeen)
                /* can't have two periods */
                flags &= ~DATA_IS_FLOAT;
              periodSeen = 1;
              if (data != data0 && !(isdigit(*(data - 1)) || *(data - 1) == '+' || *(data - 1) == '-'))
                /* period is not preceded by a number or sign */
                flags &= ~DATA_IS_FLOAT;
            }
          if ((*data == 'e' || *data == 'E') && flags & DATA_IS_FLOAT)
            {
              /* exponent if [.0123456789]e[+-0123456789] */
              if (data == data0 || (!isdigit(*(data - 1)) && *(data - 1) != '.'))
                flags &= ~DATA_IS_FLOAT;
              else if (*(data + 1) == '+' || *(data + 1) == '-' || isdigit(*(data + 1)))
                {
                  data++;
                  exponentialSeen = 1;
                }
              else
                flags &= ~DATA_IS_FLOAT;
            }
        }
      data++;
    }
  while (*data && exponentialSeen)
    {
      if (!isdigit(*data))
        {
          flags = 0;
          break;
        }
      data++;
    }
  if ((flags & DATA_IS_INT) && sscanf(data0, "%ld", &longValue) == 1)
    return SDDS_LONG;
  if ((flags & DATA_IS_FLOAT) && sscanf(data0, "%lf", &doubleValue) == 1)
    return SDDS_DOUBLE;
  return SDDS_STRING;
}

#define COLUMN_MACRO_MODE 0
#define PARAMETER_MACRO_MODE 1
#define DICTIONARY_MACRO_MODE 2
#define N_MACRO_MODES 3

void doMacroOutput(NAMELIST_TEXT *nl_text, FIELDS *EM_problem, char **macroTag, char **macroValue, long macros)
{
  char *outputFilename;
  char *modeOption[N_MACRO_MODES] = {"columns", "parameters", "dictionary"};
  SDDS_DATASET SDDSout;
  long i, *type, longValue;
  double doubleValue;

  macro_output_struct.filename = NULL;
  macro_output_struct.mode = "parameters";

  process_namelist(&macro_output, nl_text);
  print_namelist(stdout, &macro_output);

  if (!macro_output_struct.filename)
    bomb("no filename given for macro_output", NULL);

  outputFilename = compose_filename(macro_output_struct.filename, EM_problem->rootname);
  type = tmalloc(sizeof(*type) * macros);

  if (!SDDS_InitializeOutput(&SDDSout, SDDS_BINARY, 0, NULL, "spiffe macro output", outputFilename))
    {
      fprintf(stderr, "Problem creating macro output file %s\n", outputFilename);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }

  switch (match_string(macro_output_struct.mode, modeOption, N_MACRO_MODES, 0))
    {
    case COLUMN_MACRO_MODE:
      for (i = 0; i < macros; i++)
        {
          type[i] = identifyTypeOfData(macroValue[i]);
          if (!SDDS_DefineSimpleColumn(&SDDSout, macroTag[i], NULL, type[i]))
            {
              SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
              exit(1);
            }
        }
      if (!SDDS_WriteLayout(&SDDSout) || !SDDS_StartPage(&SDDSout, 1))
        {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          exit(1);
        }
      for (i = 0; i < macros; i++)
        {
          switch (type[i])
            {
            case SDDS_LONG:
              if (!sscanf(macroValue[i], "%ld", &longValue) ||
                  !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, 0, i, longValue, -1))
                {
                  fprintf(stderr, "Problem setting or scanning value %s for macro %s\n",
                          macroValue[i], macroTag[i]);
                  SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
                  exit(1);
                }
              break;
            case SDDS_DOUBLE:
              if (!sscanf(macroValue[i], "%lf", &doubleValue) ||
                  !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, 0, i, doubleValue, -1))
                {
                  fprintf(stderr, "Problem setting or scanning value %s for macro %s\n",
                          macroValue[i], macroTag[i]);
                  SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
                  exit(1);
                }
              break;
            case SDDS_STRING:
              if (!SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, 0, i, macroValue[i], -1))
                {
                  fprintf(stderr, "Problem setting value %s for macro %s\n",
                          macroValue[i], macroTag[i]);
                  SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
                  exit(1);
                }
              break;
            default:
              break;
            }
        }
      break;
    case PARAMETER_MACRO_MODE:
      for (i = 0; i < macros; i++)
        {
          type[i] = identifyTypeOfData(macroValue[i]);
          if (!SDDS_DefineSimpleParameter(&SDDSout, macroTag[i], NULL, type[i]))
            {
              SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
              exit(1);
            }
        }
      if (!SDDS_WriteLayout(&SDDSout) || !SDDS_StartPage(&SDDSout, 1))
        {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          exit(1);
        }
      for (i = 0; i < macros; i++)
        {
          switch (type[i])
            {
            case SDDS_LONG:
              if (!sscanf(macroValue[i], "%ld", &longValue) ||
                  !SDDS_SetParameters(&SDDSout, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, i, longValue, -1))
                {
                  fprintf(stderr, "Problem setting or scanning value %s for macro %s\n",
                          macroValue[i], macroTag[i]);
                  SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
                  exit(1);
                }
              break;
            case SDDS_DOUBLE:
              if (!sscanf(macroValue[i], "%lf", &doubleValue) ||
                  !SDDS_SetParameters(&SDDSout, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, i, doubleValue, -1))
                {
                  fprintf(stderr, "Problem setting or scanning value %s for macro %s\n",
                          macroValue[i], macroTag[i]);
                  SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
                  exit(1);
                }
              break;
            case SDDS_STRING:
              if (!SDDS_SetParameters(&SDDSout, SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE, i, macroValue[i], -1))
                {
                  fprintf(stderr, "Problem setting value %s for macro %s\n",
                          macroValue[i], macroTag[i]);
                  SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
                  exit(1);
                }
              break;
            default:
              fprintf(stderr, "invalid data type %ld writing macro file (%s=%s)\n",
                      type[i], macroTag[i], macroValue[i]);
              exit(1);
              break;
            }
        }
      break;
    case DICTIONARY_MACRO_MODE:
      if (!SDDS_DefineSimpleColumn(&SDDSout, "MacroTag", NULL, SDDS_STRING) ||
          !SDDS_DefineSimpleColumn(&SDDSout, "MacroStringValue", NULL, SDDS_STRING))
        {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          exit(1);
        }
      if (!SDDS_WriteLayout(&SDDSout) || !SDDS_StartPage(&SDDSout, macros) ||
          !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_INDEX, macroTag, macros, 0) ||
          !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_INDEX, macroValue, macros, 1))
        {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          exit(1);
        }
      break;
    default:
      fprintf(stderr, "unknown mode \"%s\" for macro_output\n", macro_output_struct.mode);
      exit(1);
      break;
    }

  if (!SDDS_WritePage(&SDDSout) || !SDDS_Terminate(&SDDSout))
    {
      fprintf(stderr, "Problem writing macro_output file %s\n", outputFilename);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }

  free(type);
}
