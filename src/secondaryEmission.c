/*************************************************************************\
 * Copyright (c) 2005 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2005 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: secondary_emission.c
 * purpose: setup and execution of secondary emission
 *
 * Michael Borland, 2005
 */
#include "spiffe.h"

#define DEBUG 0

long inversePoissonCDF(double mu, double C);

void processSecondaryEmissionDefinition(
                                        SECONDARY_EMISSION **secondaryEmission,
                                        long *nSecondaryDefns,
                                        NAMELIST_TEXT *nl_text,
                                        FIELDS *EM)
{
#include "secondaryEmission.h"
  long i;
  SDDS_DATASET SDDSinput;
  SECONDARY_EMISSION *SE;
  double energyMin, energyMax;

  input_file = kinetic_energy_column = yield_column = NULL;
  emitted_momentum = 0;
  yield_limit = 0;
  verbosity = 1;
  material_id = 1;
  log_file = NULL;

  process_namelist(&secondary_emission, nl_text);
  print_namelist(stdout, &secondary_emission);

  if (!input_file)
    bomb("supply an input file", NULL);
  if (!kinetic_energy_column)
    bomb("supply the name of a column containing energy values", NULL);
  if (!yield_column)
    bomb("supply the name of a column containing yield values", NULL);
  if (emitted_momentum < 0)
    bomb("emitted_momentum < 0", NULL);

  if (!(*secondaryEmission = SDDS_Realloc(*secondaryEmission, sizeof(**secondaryEmission) * (*nSecondaryDefns + 1))))
    bomb("memory allocation failure", NULL);

  SE = *secondaryEmission + *nSecondaryDefns;

  SE->emittedMomentum = emitted_momentum;
  SE->verbosity = verbosity;
  SE->materialID = material_id;
  SE->yieldLimit = yield_limit;
  SE->generation = SE->directionCode = NULL;
  SE->z = SE->r = SE->p = SE->Q = SE->yieldData = SE->energyData = NULL;
  SE->particleID = NULL;
  SE->logFile = NULL;
  SE->maxQueueLength = SE->queueLength = 0;
  SE->nData = 0;
  if (!SDDS_CopyString(&SE->inputFile, input_file))
    SDDS_Bomb("Memory allocation failure");
  if (!SDDS_InitializeInput(&SDDSinput, input_file) ||
      SDDS_ReadPage(&SDDSinput) != 1)
    SDDS_Bomb("problem reading secondary emission yield data file");
  if ((SE->nData = SDDS_CountRowsOfInterest(&SDDSinput)) <= 0)
    bomb("no data in yield file", NULL);
  if (SDDS_CheckColumn(&SDDSinput, kinetic_energy_column, "eV", SDDS_ANY_FLOATING_TYPE,
                       stderr) != SDDS_CHECK_OKAY)
    SDDS_Bomb("Problem with kinetic_energy_column (check existence, type, and units (eV)).");
  if (SDDS_CheckColumn(&SDDSinput, yield_column, "", SDDS_ANY_FLOATING_TYPE,
                       stderr) != SDDS_CHECK_OKAY)
    SDDS_Bomb("Problem with yield_column (check existence and type");

  if (!(SE->energyData = SDDS_GetColumnInDoubles(&SDDSinput, kinetic_energy_column)) ||
      !(SE->yieldData = SDDS_GetColumnInDoubles(&SDDSinput, yield_column)))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors | SDDS_VERBOSE_PrintErrors);
  find_min_max(&energyMin, &energyMax, SE->energyData, SE->nData);
  SE->pMin = sqrt(sqr(energyMin / me_mev / 1e6 + 1) - 1);
  SE->pMax = sqrt(sqr(energyMax / me_mev / 1e6 + 1) - 1);

  for (i = 1; i < SE->nData; i++)
    if (SE->energyData[i - 1] >= SE->energyData[i])
      bomb("secondary emission data is not ordered by increasing energy values", NULL);

  if (log_file)
    {
      char label[16384], *filename;
      SDDS_CopyString(&SE->logFile, log_file);
      filename = compose_filename(log_file, EM->rootname);
      sprintf(label, "Primary/secondary activity due to material %s\n", SE->inputFile);
      if (!SDDS_InitializeOutput(&SE->SDDSlog, SDDS_BINARY, 0, label, NULL, filename) ||
          !SDDS_DefineSimpleColumn(&SE->SDDSlog, "z", "m", SDDS_FLOAT) ||
          !SDDS_DefineSimpleColumn(&SE->SDDSlog, "r", "m", SDDS_FLOAT) ||
          !SDDS_DefineSimpleColumn(&SE->SDDSlog, "t", "s", SDDS_FLOAT) ||
          !SDDS_DefineSimpleColumn(&SE->SDDSlog, "K", "eV", SDDS_FLOAT) ||
          !SDDS_DefineSimpleColumn(&SE->SDDSlog, "yield", NULL, SDDS_FLOAT) ||
          !SDDS_DefineSimpleColumn(&SE->SDDSlog, "emitted", NULL, SDDS_SHORT) ||
          !SDDS_DefineSimpleColumn(&SE->SDDSlog, "directionCode", NULL, SDDS_SHORT) ||
          !SDDS_DefineSimpleColumn(&SE->SDDSlog, "generation", NULL, SDDS_SHORT) ||
          !SDDS_DefineSimpleColumn(&SE->SDDSlog, "particleID", NULL, SDDS_LONG) ||
          !SDDS_WriteLayout(&SE->SDDSlog))
        {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          exit(1);
        }
      free(filename);
    }

  SE->initialized = 1;
  *nSecondaryDefns += 1;
}

void queueSecondaryEmission(SECONDARY_EMISSION *secondaryEmission, long nSecondaryDefns,
                            double p, double z, double r, short directionCode, double Q,
                            short generation, FIELDS *EM_problem, BEAM *beam)
{
  long iz, ir, iz1, ir1, im;
  short material;
  SECONDARY_EMISSION *SE;
  double Ez, Er, KSecondary;

  if ((z + EM_problem->dz) >= EM_problem->zmax)
    return;

  iz = (z - EM_problem->zmin) / EM_problem->dz;
  ir = fabs(r) / EM_problem->dr;
  iz1 = iz + 1;
  ir1 = ir + 1;
  if (iz >= EM_problem->nz || ir >= EM_problem->nr || iz < 0 || ir < 0)
    return;
  if (EM_problem->metal_flags[iz][ir] & FL_IS_INSIDE)
    return;
  else if (EM_problem->metal_flags[iz][ir] & FL_IS_SURFACE)
    {
      if (iz1 > EM_problem->nz)
        return;
      else
        {
          if (ir1 > EM_problem->nr)
            return;
          else
            {
              if (EM_problem->metal_flags[iz1][ir1] & FL_IS_METAL &&
                  EM_problem->metal_flags[iz][ir1] & FL_IS_METAL &&
                  EM_problem->metal_flags[iz1][ir] & FL_IS_METAL)
                return;
            }
        }
    }

  if ((material = EM_problem->materialID[iz][ir]) < 1)
    {
      if ((material = EM_problem->materialID[iz1][ir]) < 1)
        {
          if ((material = EM_problem->materialID[iz][ir1]) < 1)
            {
              if ((material = EM_problem->materialID[iz1][ir1]) < 1)
                {
                  /* this location has no secondary emission */
                  return;
                }
            }
        }
    }

  for (im = 0; im < nSecondaryDefns; im++)
    {
      if (material == secondaryEmission[im].materialID)
        break;
    }
  if (im == nSecondaryDefns)
    {
      fprintf(stderr, "Error: no match for material %hd  for r=%e, z=%e, (ir=%ld, iz=%ld)\n",
              material, r, z, ir, iz);
      exit(1);
    }
  SE = secondaryEmission + im;

  if (SE->pMin > p || SE->pMax < p)
    return;

  interpolateFields(&Ez, &Er, NULL, NULL, NULL, NULL, r, z, EM_problem);
  KSecondary = me_mev * 1e6 * (sqrt(sqr(SE->emittedMomentum) + 1) - 1);
  if ((directionCode == GOING_LEFT && (Ez * EM_problem->dz * 1e-9) > KSecondary) ||
      (directionCode == GOING_RIGHT && (-1 * Ez * EM_problem->dz * 1e-9) > KSecondary) ||
      (directionCode == GOING_DOWN && (Er * EM_problem->dr * 1e-9) > KSecondary) ||
      (directionCode == GOING_UP && (-1 * Er * EM_problem->dr * 1e-9) > KSecondary))
    return;

  if (SE->queueLength >= SE->maxQueueLength)
    {
      SE->maxQueueLength += 100;
      if (SE->verbosity > 0)
        {
          printf("Increasing size of secondary emission queue for material %hd to %ld\n",
                 material, SE->maxQueueLength);
          fflush(stdout);
        }
      if (!(SE->z = SDDS_Realloc(SE->z, sizeof(*(SE->z)) * SE->maxQueueLength)) ||
          !(SE->r = SDDS_Realloc(SE->r, sizeof(*(SE->r)) * SE->maxQueueLength)) ||
          !(SE->p = SDDS_Realloc(SE->p, sizeof(*(SE->p)) * SE->maxQueueLength)) ||
          !(SE->Q = SDDS_Realloc(SE->Q, sizeof(*(SE->Q)) * SE->maxQueueLength)) ||
          !(SE->directionCode =
            SDDS_Realloc(SE->directionCode,
                         sizeof(*(SE->directionCode)) * SE->maxQueueLength)) ||
          !(SE->generation =
            SDDS_Realloc(SE->generation,
                         sizeof(*(SE->generation)) * SE->maxQueueLength)) ||
          !(SE->particleID =
            SDDS_Realloc(SE->particleID,
                         sizeof(*(SE->particleID)) * SE->maxQueueLength)))
        SDDS_Bomb("memory allocation failure increasing size of secondary emission queue");
    }

  if (SE->verbosity > 1)
    {
      printf("Adding particles to secondary emission queue for material %hd at z=%e, r=%e p=%e, Q=%e, gen=%hd\n",
             material, z, r, p, Q, generation);
      fflush(stdout);
    }

  SE->z[SE->queueLength] = z;
  SE->r[SE->queueLength] = r;
  SE->p[SE->queueLength] = p;
  SE->Q[SE->queueLength] = Q;
  SE->directionCode[SE->queueLength] = directionCode;
  SE->generation[SE->queueLength] = generation;
  SE->particleID[SE->queueLength] = beam->nextParticleID++;
  SE->queueLength += 1;
}

void findLossCoordinates(double *zReturn, double *rReturn, short *directionCode,
                         double z0, double r0, double z1, double r1,
                         FIELDS *EM_problem)
{
  double zm, rm;
  long iz1, ir1, iz0, ir0;

  *directionCode = 0;
  iz0 = (z0 - EM_problem->zmin) / EM_problem->dz;
  ir0 = r0 / EM_problem->dr;
  if (z1 < EM_problem->zmin)
    iz1 = -1;
  else
    iz1 = (z1 - EM_problem->zmin) / EM_problem->dz;
  if (r1 < 0)
    ir1 = -1;
  else
    ir1 = r1 / EM_problem->dr;

  if (r0 == r1 && z0 == z1)
    {
      *directionCode = 1;
      *zReturn = z0;
      *rReturn = r0;
      return;
    }

  if (ir0 == ir1)
    {
      /* hit a constant z surface */
      if (iz0 > iz1)
        *directionCode = GOING_LEFT;
      else
        *directionCode = GOING_RIGHT;
    }
  else if (iz0 == iz1)
    {
      /* hit a constant r surface */
      if (ir0 > ir1)
        *directionCode = GOING_DOWN;
      else
        *directionCode = GOING_UP;
    }
  else
    {
      if (z1 < z0)
        {
          /* traveling left */
          if (z1 < 0)
            *directionCode = GOING_LEFT; /* hit a constant z surface */
          else if (ir0 != ir1)
            {
              if (ir1 >= 0 && ir1 < EM_problem->nr && EM_problem->metal_flags[iz0][ir1] & FL_IS_METAL)
                *directionCode = GOING_LEFT;
              else if (ir1 < 0 && EM_problem->metal_flags[iz0][-ir1] & FL_IS_METAL)
                *directionCode = GOING_LEFT;
            }
          else if (ir0 == ir1 && (ir0 < (EM_problem->nr - 1) && EM_problem->metal_flags[iz0][ir0 + 1] & FL_IS_METAL))
            *directionCode = GOING_LEFT;
          else
            {
              /* assume it hit a constant r surface */
              if (r1 < r0)
                *directionCode = GOING_DOWN;
              else
                *directionCode = GOING_UP;
            }
        }
      else if (z1 > z0)
        {
          /* traveling left */
          if (z1 > EM_problem->zmax)
            *directionCode = GOING_RIGHT; /* hit a constant z surface */
          else if (ir0 != ir1)
            {
              if (ir1 >= 0 && ir1 < EM_problem->nr && EM_problem->metal_flags[iz1][ir1] & FL_IS_METAL)
                *directionCode = GOING_RIGHT;
              else if (ir1 < 0 && EM_problem->metal_flags[iz1][-ir1] & FL_IS_METAL)
                *directionCode = GOING_RIGHT;
            }
          else if (ir0 == ir1 && (ir0 < (EM_problem->nr - 1) && EM_problem->metal_flags[iz1][ir0 + 1] & FL_IS_METAL))
            *directionCode = GOING_RIGHT;
          else
            {
              /* assume it hit a constant r surface */
              if (r1 < r0)
                *directionCode = GOING_DOWN;
              else
                *directionCode = GOING_UP;
            }
        }
      else
        {
          if (r1 < r0)
            *directionCode = GOING_DOWN;
          else
            *directionCode = GOING_UP;
        }
    }

  if (*directionCode == GOING_LEFT)
    {
      zm = ((long)((z0 - EM_problem->zmin) / EM_problem->dz)) * EM_problem->dz + EM_problem->zmin + EM_problem->dz * 1e-12;
      rm = (r1 - r0) / (z1 - z0) * (zm - z0) + r0;
    }
  else if (*directionCode == GOING_RIGHT)
    {
      zm = ((long)((z0 - EM_problem->zmin) / EM_problem->dz)) * EM_problem->dz + EM_problem->zmin + EM_problem->dz * (1 - 1e-12);
      rm = (r1 - r0) / (z1 - z0) * (zm - z0) + r0;
    }
  else if (*directionCode == GOING_DOWN)
    {
      rm = ((long)(r0 / EM_problem->dr)) * EM_problem->dr + EM_problem->dr * 1e-12;
      zm = (z1 - z0) / (r1 - r0) * (rm - r0) + z0;
    }
  else
    {
      rm = ((long)(r0 / EM_problem->dr)) * EM_problem->dr + EM_problem->dr * (1 - 1e-12);
      zm = (z1 - z0) / (r1 - r0) * (rm - r0) + z0;
    }

  /*
    fprintf(stdout, "z0=%e, r0=%e, z1=%e, r1=%e, iz0=%ld, iz1=%ld, ir0=%ld, ir1=%ld, code=%c, zm=%e, rm=%e\n",
    z0, r0, z1, r1, iz0, iz1, ir0, ir1, 
    (*directionCode==GOING_DOWN?'d':
    (*directionCode==GOING_UP?'u':
    (*directionCode==GOING_LEFT?'l':
    (*directionCode==GOING_RIGHT?'r':'?')))),
    zm, rm);
  */

  *zReturn = zm;
  *rReturn = fabs(rm);
}

void emitSecondaryElectrons(BEAM *beam, SECONDARY_EMISSION *secondaryEmission, long nSecondaryDefns,
                            FIELDS *EM_problem)
{
  long iSource, im;
  long ip, np0, npNew, interpCode;
  double meanYield, gamma, kineticEnergy;
  double z, r, theta;
  long total;
  SECONDARY_EMISSION *SE;
  static float *zLog = NULL, *rLog = NULL, *tLog = NULL, *KLog = NULL, *YieldLog = NULL;
  static short *EmittedLog = NULL, *DirectionLog = NULL, *generationLog = NULL;
  static int32_t *particleIDLog = NULL;
  static long logBufferSize = 0;
  long nLogged;

  for (im = 0; im < nSecondaryDefns; im++)
    {
      SE = secondaryEmission + im;
      nLogged = 0;

      if (!SE->initialized || !SE->queueLength)
        continue;

      if (SE->verbosity > 0)
        {
          printf("%ld secondary emission sources for material %hd\n", SE->queueLength,
                 SE->materialID);
          fflush(stdout);
        }

      np0 = beam->np;
      for (iSource = 0; iSource < SE->queueLength; iSource++)
        {
          r = SE->r[iSource];
          z = SE->z[iSource];
          if (SE->verbosity > 2)
            {
              printf("Secondary location for source %ld: z=%e, r=%e\n",
                     iSource, z, r);
              printf("Source properties: p=%e, Q=%e\n", SE->p[iSource],
                     SE->Q[iSource]);
              fflush(stdout);
            }

          /* Compute number of particles to be emitted */
          gamma = sqrt(sqr(SE->p[iSource]) + 1);
          kineticEnergy = (gamma - 1) * me_mev * 1e6;
          meanYield = interp(SE->yieldData, SE->energyData,
                             SE->nData, kineticEnergy, 0, 1, &interpCode);
          if (!interpCode)
            {
              if (SE->verbosity > 0)
                {
                  printf("Secondary yield curve interpolation failed for K=%e eV\n", kineticEnergy);
                  fflush(stdout);
                }
              continue;
            }
          if (SE->verbosity > 2)
            {
              printf("Secondary yield K=%e (p=%e, gamma=%e) is %e\n", kineticEnergy,
                     SE->p[iSource], gamma, meanYield);
              fflush(stdout);
            }

          /* Use poisson distribution to determine number of particles to emit */
          total = inversePoissonCDF(meanYield, random_1(1));
          if (SE->yieldLimit > 0 && total > SE->yieldLimit)
            total = SE->yieldLimit;

          /* Emit particles if needed */
          if (total > 0)
            {
              /* Resize particle arrays if needed */
              npNew = beam->np + total;
              if (npNew > beam->max_np)
                {
                  beam->max_np = 1.1 * (npNew + 1);
                  if (beam->max_np > LONG_MAX / 2 || beam->max_np <= 0)
                    bomb("integer overflow in computing particle array dimension", NULL);
                  beam->z = trealloc(beam->z, sizeof(*beam->z) * beam->max_np);
                  beam->r = trealloc(beam->r, sizeof(*beam->r) * beam->max_np);
                  beam->pz = trealloc(beam->pz, sizeof(*beam->pz) * beam->max_np);
                  beam->pr = trealloc(beam->pr, sizeof(*beam->pr) * beam->max_np);
                  beam->vz = trealloc(beam->vz, sizeof(*beam->vz) * beam->max_np);
                  beam->vr = trealloc(beam->vr, sizeof(*beam->vr) * beam->max_np);
                  beam->r0 = trealloc(beam->r0, sizeof(*beam->r0) * beam->max_np);
                  beam->z0 = trealloc(beam->z0, sizeof(*beam->z0) * beam->max_np);
                  beam->t0 = trealloc(beam->t0, sizeof(*beam->t0) * beam->max_np);
                  beam->rmin = trealloc(beam->rmin, sizeof(*beam->rmin) * beam->max_np);
                  beam->zmin = trealloc(beam->zmin, sizeof(*beam->zmin) * beam->max_np);
                  beam->rmax = trealloc(beam->rmax, sizeof(*beam->rmax) * beam->max_np);
                  beam->zmax = trealloc(beam->zmax, sizeof(*beam->zmax) * beam->max_np);
                  beam->r_pphi = trealloc(beam->r_pphi, sizeof(*beam->r_pphi) * beam->max_np);
                  beam->vphi = trealloc(beam->vphi, sizeof(*beam->vphi) * beam->max_np);
                  beam->Q = trealloc(beam->Q, sizeof(*beam->Q) * beam->max_np);
                  beam->generation = trealloc(beam->generation, sizeof(*beam->generation) * beam->max_np);
                  beam->particleID = trealloc(beam->particleID, sizeof(*beam->particleID) * beam->max_np);
                  beam->status = trealloc(beam->status, sizeof(*beam->status) * beam->max_np);
                }

              for (ip = beam->np; ip < npNew; ip++)
                {
                  beam->generation[ip] = SE->generation[iSource] + 1;
                  beam->particleID[ip] = SE->particleID[iSource];
                  beam->Q[ip] = SE->Q[iSource];
                  switch (SE->directionCode[iSource])
                    {
                    case GOING_RIGHT: /* secondary goes left */
                      theta = -random_1(1) * PI;
                      break;
                    case GOING_LEFT: /* secondary goes right */
                      theta = random_1(1) * PI;
                      break;
                    case GOING_UP: /* secondary goes down */
                      theta = random_1(1) * PI + PI / 2;
                      break;
                    case GOING_DOWN: /* secondary goes up */
                      theta = random_1(1) * PI - PI / 2;
                      break;
                    default:
                      theta = random_1(1) * PIx2;
                      break;
                    }
                  beam->rmin[ip] = beam->zmin[ip] = DBL_MAX;
                  beam->rmax[ip] = beam->zmax[ip] = -DBL_MAX;
                  beam->pz[ip] = SE->emittedMomentum * sin(theta);
                  beam->pr[ip] = SE->emittedMomentum * cos(theta);
                  gamma = sqrt(sqr(beam->pz[ip]) + sqr(beam->pr[ip]) + 1);
                  beam->vz[ip] = c_mks * beam->pz[ip] / gamma;
                  beam->vr[ip] = c_mks * beam->pr[ip] / gamma;
                  beam->z[ip] = beam->z0[ip] = z;
                  beam->r[ip] = beam->r0[ip] = r;
                  beam->t0[ip] = EM_problem->time;
                  beam->vphi[ip] = beam->r_pphi[ip] = 0;
                  beam->status[ip] = PART_ACTIVE;
                  beam->npActive++;
                }

              if (SE->verbosity > 2)
                {
                  printf("%ld particles present after secondary source %ld\n",
                         npNew, iSource);
                  fflush(stdout);
                }

              if (SE->logFile)
                {
                  if (nLogged >= logBufferSize)
                    {
                      logBufferSize += 10;
                      if (!(zLog = SDDS_Realloc(zLog, sizeof(*zLog) * logBufferSize)) ||
                          !(rLog = SDDS_Realloc(rLog, sizeof(*rLog) * logBufferSize)) ||
                          !(tLog = SDDS_Realloc(tLog, sizeof(*tLog) * logBufferSize)) ||
                          !(KLog = SDDS_Realloc(KLog, sizeof(*KLog) * logBufferSize)) ||
                          !(YieldLog = SDDS_Realloc(YieldLog, sizeof(*YieldLog) * logBufferSize)) ||
                          !(generationLog = SDDS_Realloc(generationLog, sizeof(*generationLog) * logBufferSize)) ||
                          !(particleIDLog = SDDS_Realloc(particleIDLog, sizeof(*particleIDLog) * logBufferSize)) ||
                          !(EmittedLog = SDDS_Realloc(EmittedLog, sizeof(*EmittedLog) * logBufferSize)) ||
                          !(DirectionLog = SDDS_Realloc(DirectionLog, sizeof(*DirectionLog) * logBufferSize)))
                        SDDS_Bomb("memory allocation failure (secondary logging)");
                    }
                  zLog[nLogged] = z;
                  rLog[nLogged] = r;
                  tLog[nLogged] = EM_problem->time;
                  KLog[nLogged] = kineticEnergy;
                  YieldLog[nLogged] = meanYield;
                  EmittedLog[nLogged] = npNew - beam->np;
                  DirectionLog[nLogged] = SE->directionCode[iSource];
                  generationLog[nLogged] = SE->generation[iSource] + 1;
                  particleIDLog[nLogged] = beam->nextParticleID++;
                  nLogged++;
                }

              beam->np = npNew;
            }
        }

      if (SE->logFile && nLogged)
        {
          if (!SDDS_StartPage(&SE->SDDSlog, nLogged) ||
              !SDDS_SetColumn(&SE->SDDSlog, SDDS_SET_BY_NAME, zLog, nLogged, "z") ||
              !SDDS_SetColumn(&SE->SDDSlog, SDDS_SET_BY_NAME, rLog, nLogged, "r") ||
              !SDDS_SetColumn(&SE->SDDSlog, SDDS_SET_BY_NAME, tLog, nLogged, "t") ||
              !SDDS_SetColumn(&SE->SDDSlog, SDDS_SET_BY_NAME, KLog, nLogged, "K") ||
              !SDDS_SetColumn(&SE->SDDSlog, SDDS_SET_BY_NAME, YieldLog, nLogged, "yield") ||
              !SDDS_SetColumn(&SE->SDDSlog, SDDS_SET_BY_NAME, EmittedLog, nLogged, "emitted") ||
              !SDDS_SetColumn(&SE->SDDSlog, SDDS_SET_BY_NAME, DirectionLog, nLogged, "directionCode") ||
              !SDDS_SetColumn(&SE->SDDSlog, SDDS_SET_BY_NAME, generationLog, nLogged, "generation") ||
              !SDDS_SetColumn(&SE->SDDSlog, SDDS_SET_BY_NAME, generationLog, nLogged, "particleID") ||
              !SDDS_WritePage(&SE->SDDSlog))
            {
              SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
              exit(1);
            }
        }
      nLogged = 0;

      if (SE->verbosity > 0)
        {
          printf("Secondary emission increased number of particles from %ld to %ld\n",
                 np0, beam->np);
          fflush(stdout);
        }
      SE->queueLength = 0;
    }
}

long inversePoissonCDF(double mu, double C)
{
  double sum, expMinusMu, term;
  long r, rMax;

  r = 0;
  rMax = 10 * mu;
  expMinusMu = exp(-mu);
  sum = expMinusMu;
  term = 1;
  while (r <= rMax && C > sum)
    {
      term *= mu / (++r);
      sum += term * expMinusMu;
    }
  return r;
}
