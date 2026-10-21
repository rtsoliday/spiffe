/* file: semaphores.c
 * purpose: creation of semaphore files indicating that a run has started,
 *          completed, or failed.  Modeled on the equivalent feature in elegant.
 */
#include "spiffe.h"
#include "semaphores.h"

static char *doneSemaphoreFile = NULL;
static char *failedSemaphoreFile = NULL;
static short doneSemaphoreCreated = 0;

void createSemaphoreFile(char *filename)
{
  FILE *fp;

  if (!filename || !strlen(filename) || is_blank(filename))
    return;
  printf("Creating semaphore file \"%s\"\n", filename);
  fflush(stdout);
  if (!(fp = fopen(filename, "w")))
    {
      fprintf(stderr, "Problem creating semaphore file \"%s\"\n", filename);
      exit(1);
    }
  if (wild_match(filename, "*.done"))
    fprintf(fp, "%8.2f\n", cpu_time() / 100.0);
  fclose(fp);
}

void createDoneSemaphoreFile(void)
{
  doneSemaphoreCreated = 1;
  createSemaphoreFile(doneSemaphoreFile);
}

static void createFailedSemaphoreFile(void)
{
  /* Registered with atexit(), so it runs on normal return from main() as well
   * as on every exit() and bomb() path.
   */
  static short inHandler = 0;
  if (inHandler || doneSemaphoreCreated)
    return;
  inHandler = 1;
  createSemaphoreFile(failedSemaphoreFile);
}

static char *composeSemaphoreName(char *template, char *rootname)
{
  char *writable, *name, *result;

  /* compose_filename temporarily writes into the template, so it must not be
   * given the string literal that serves as the namelist default.
   */
  SDDS_CopyString(&writable, template);
  name = compose_filename(writable, rootname);
  SDDS_CopyString(&result, name);
  if (name != writable)
    free(name);
  free(writable);
  return result;
}

void process_semaphores(NAMELIST_TEXT *nl_text, FIELDS *EM_problem)
{
  static short handlerRegistered = 0;
  char *startedSemaphoreFile;

  started = "%s.started";
  done = "%s.done";
  failed = "%s.failed";
  immediate = NULL;

  process_namelist(&semaphores, nl_text);
  print_namelist(stdout, &semaphores);

  if (immediate && strlen(immediate))
    {
      immediate = composeSemaphoreName(immediate, EM_problem->rootname);
      createSemaphoreFile(immediate);
      free(immediate);
      immediate = NULL;
    }

  if (doneSemaphoreFile)
    free(doneSemaphoreFile);
  if (failedSemaphoreFile)
    free(failedSemaphoreFile);
  doneSemaphoreFile = failedSemaphoreFile = NULL;

  if (done && strlen(done))
    doneSemaphoreFile = composeSemaphoreName(done, EM_problem->rootname);
  if (failed && strlen(failed))
    failedSemaphoreFile = composeSemaphoreName(failed, EM_problem->rootname);

  if (doneSemaphoreFile && fexists(doneSemaphoreFile))
    remove(doneSemaphoreFile);
  if (failedSemaphoreFile && fexists(failedSemaphoreFile))
    remove(failedSemaphoreFile);

  if (started && strlen(started))
    {
      startedSemaphoreFile = composeSemaphoreName(started, EM_problem->rootname);
      if (fexists(startedSemaphoreFile))
        remove(startedSemaphoreFile);
      createSemaphoreFile(startedSemaphoreFile);
      free(startedSemaphoreFile);
    }

  if (failedSemaphoreFile && !handlerRegistered)
    {
      atexit(createFailedSemaphoreFile);
      handlerRegistered = 1;
    }
}
