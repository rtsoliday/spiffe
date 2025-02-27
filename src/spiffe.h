/* file: spiffe.h
 * purpose: main include file for particle-in-cell code spiffe
 *
 * Michael Borland, 1992 
 */
#include "mdb.h"
#include "namelist.h"
#include "SDDS.h"

/* 
 * Structure ON_AXIS_FIELD holds information specifying a field expansion from
 * on-axis field data.
 * See M. Borland's thesis, 2.2.7.
 */

typedef struct {
  long points, expansion_order;
  double phase, frequency;
  double *z, *Ez, *DzEz, *Dz2Ez, *Dz3Ez;
} ON_AXIS_FIELD;

/*
 * Structure FIELDS holds complete information on a solution of Maxwell's equations 
 * at a given time, plus external fields.
 *
 */

typedef struct {
#define FORMAT_VERSION -3
    long version;             /* used for field saves/readbacks */
    /* definition of spatial grid: */
    int32_t nz, nr;              /* numbers of metal grid points in z, r */
    double zmin, zmax, rmax;
    double dz, dr;            /* dX = (Xmax - Xmin)/(nX - 1) */

    /* boundary conditions */
    long lower_bc, upper_bc, right_bc, left_bc;
#define DIRICHLET 0
#define NEUMANN   1

    /* definition of time grid: */
    double start_time;        /* time of solution start */
    double time;              /* current value of time */
    double dti;               /* time step for integration */
    long time_step;           /* number of time steps taken during this run */
    long n_time_steps;        /* total number of steps that will be taken */
    long status_interval;     /* interval in steps at which to print status information */
    
    unsigned long modes;               /* flag word for field "modes", other modes */
#define FL_TM_FIELDS 1
#define FL_TE_FIELDS 2
#define FL_SPACE_CHARGE 4
#define FL_RADIAL_INTERP 8
#define FL_LONGIT_INTERP 16
#define FL_RADIAL_SMEAR 32
#define FL_LONGIT_SMEAR 64
    double zt;                /* amount by which problem has been translated */

    /* TM mode fields */ 
    double **Ez, **Jz;    /* Ez[iz][ir]   is value at z=zmin+(iz-1/2)*dz, r=ir*dr       , t=it*dt      , 0<=iz<nz+1, 0<=ir<nr    */
    double **Er, **Jr;    /* Er[iz][ir]   is value at z=zmin+iz*dz      , r=(ir-1/2)*dr , t=it*dt      , 0<=iz<nz  , 0<=ir<nr+1  */
    double **Bphi;        /* Bphi[iz][ir] is value at z=zmin+(iz-1/2)*dz, r=(ir-1/2)*dr , t=(it+1/2)*dt, 0<=iz<nz+1, 0<=ir<nr+1  */
    double **Bphi2;       /* Bphi2 is Bphi at t=it*dt, i.e., at same time as E fields */

    /* TE mode fields */
    double **Ephi, **Jphi; 
                         /* Ephi[iz][ir] is value at z=zmin+iz*dz      , r=ir*dr      , t = it*dt      ,  0<=iz<nz  , 0<=ir<nr   */
    double **Br;         /* Br[iz][ir]   is value at z=zmin+(iz-1/2)*dz, r=ir*dr      , t = (it+1/2)*dt,  0<=iz<nz+1, 0<=ir<nr   */
    double **Br2;        /* Br at t=it*dt */
    double **Bz;         /* Bz[iz][ir]   is value at z=zmin+iz*dz      , r=(ir-1/2)*dr, t = (it+1/2)*dt,  0<=iz<nz  , 0<=ir<nr+1 */
    double **Bz2;        /* Bz at t=it*dt */

    /* turn off effects of fields on beam: */
    short turnOffEz, turnOffEr, turnOffEphi;
    short turnOffBz, turnOffBr, turnOffBphi;
    
    /* metal boundary/boundary condition information: */
    short **metal_flags;      /* (nz+1)x(nr+1) array to store effect of metal boundary:
                               * If metal_flags[iz][ir] and metal_flags[iz+1][ir] both have FL_IS_METAL on,
                               *     then there is a metal surface from z=zmin+iz*dz to z=zmin+(iz+1)*dz at r=ir*dr.  
                               * If metal_flags[iz][ir] and metal_flags[iz][ir+1] both have FL_IS_METAL on, 
                               *     then there is a metal surface from r=ir*dr to r=(ir+1)*dr at z=iz*dz+zmin
                               * In addition:
                               *     metal_flags[iz][ir]&FL_EZ_ZERO    -->  Ez[iz][ir] is held to zero
                               *     metal_flags[iz][ir]&FL_ER_ZERO    -->  Er[iz][ir] is held to zero
                               *     metal_flags[iz][ir]&FL_BPHI_ZERO  -->  Bphi[iz][ir] is held to zero
                               * Note that a given metal_flag element stores information about 3 different, adjacent grid points on
                               * the different meshes.  This saves effort in applying the boundary conditions.
                               */ 
#define FL_IS_SURFACE  1
#define FL_IS_INSIDE   2
#define FL_EZ_ZERO     4
#define FL_ER_ZERO     8
#define FL_IS_VACUUM  16
#define FL_BPHI_ZERO  32
#define FL_EPHI_ZERO  64
#define FL_BZ_ZERO   128
#define FL_BR_ZERO   256
#define FL_IS_METAL (FL_IS_SURFACE|FL_IS_INSIDE)
#define FL_LAST_ACTIVE_PT 512

    /* variables after this location are NOT part of the field save file */
    char *run_name, *rootname;

    /* material ID number of surface points */
    short **materialID;

    /* imposed DC electric fields */
    double **EzImposed, **ErImposed;
    double **PsiImposed;
    double imposedEFieldRampTime, imposedEFieldFlatTopTime;
    
    /* imposed DC magnetic fields (solenoids) */
    /* these arrays are both such that A[iz][ir] is the value at z=iz*dz+zmin and r=ir*dr */
    double **BzImposed, **BrImposed;
    
    /* constant imposed fields */
    double constantEz, constantBz;
    double constantEr, constantBr;
    double constantEphi, constantBphi;
    
    /* on-axis fields */
    ON_AXIS_FIELD *onAxisField;
    double EzOOEMin, EzOOEMax;
    long onAxisFields;

    double **Psi;             /* Psi[iz][ir] is the value of a potential at z=zmin+iz*dz and r=ir*dr */
    double **Q;               /* Q[iz][ir] is the amount of charge belonging to the cell around z=zmin+iz*dz and r=ir*dr */
    double *area_c;           /* area_c[ir] is the area of the annulus from dr*(ir-1/2) to dr*(ir+1/2) (Centered spatial grid) */
    double *area_o;           /* area_o[ir] is the area of the annulus from dr*(ir-1) to dr*ir (Offset from spatial grid) */
    } FIELDS;



/*
 * Structure ANTENNA holds information on a driving antenna, corresponding to a
 * single define_antenna namelist:
 *
 * antenna is driven with I=current*W(t-time_offset)*sin(2*pi*frequency*t+phase)
 *
 */

typedef struct {
    double start, end;       /* starting/ending coorinate */
    long istart, iend;       /* starting/ending indices appropriate for Jz or Jr grids */
    double position;
    long iposition;

    long direction;          /* 0 (1) is z (r) direction */
    
    double current;          /* Amperes */
    double frequency;        /* Hertz */
    double phase;            /* radians */

    long n_pts;              /* number of points in waveform */
    long last_index;
    double *t_wf, *W_wf;     /* (time(seconds), amplitude) pairs for waveform */
    double time_offset;      /* seconds */
    } ANTENNA;

/*
 * Structure RESISTOR holds information on a damping resistor, defined by its
 * length, position, and conductivity.
 *
 *
 */

typedef struct {
    long n_resistors;
    double *start, *end;       /* starting/ending coorinate */
    long *istart, *iend;       /* starting/ending indices appropriate for Jz or Jr grids */
    double *position;
    long *iposition;
    long *direction;          /* 0 (1) is z (r) direction */
    double *conductivity;     /* 1/(m*ohm) */    
    } RESISTOR;

/*
 * Structure SAMPLE_DEFN holds information defining desired field samples 
 * corresponding to a single define_field_sampling namelist.
 *
 */

typedef struct {
    char *filename;
    long component_code, direction_code;
    double min_coord, max_coord, position;
    double time_interval, start_time;
    long step_interval;      /* time_interval/dt */
    long start_step;         /* time_start/dt */
    long sample_index;
    long time_sequence, fileInitialized;
    SDDS_TABLE outTable;
    } SAMPLE_DEFN;

/*
 * Structure BEAM holds complete information on an ensemble of particles
 *
 */

typedef struct {
    long np;                    /* number of macroparticles */
    long npActive;              /* number of active macroparticles */
    long max_np;                /* maximum number that arrays will currently accomodate */
    double *Q;                  /* charge = -N*e */
    double QoMC;                /* charge/(mass*c) = -e/(me*c) */
    double stiffness;           /* mass/me */
    double *z, *r;              /* positions in meters */
    double *pz, *pr;            /* components beta*gamma */
    double *r_pphi;             /* r*beta_phi*gamma */
    double *vz, *vr, *vphi;     /* velocities in m/s */
    double *r0, *z0, *t0;       /* position and time at emission */
    double *rmax, *rmin;        /* maximum and minimum r position at any time */
    double *zmax, *zmin;        /* maximum and minimum z position at any time */
    int32_t *particleID;           /* particle ID number */
    int32_t nextParticleID;        /* particle ID for next-created particle */
    short *generation;          /* zero for primary, increments by 1 for each secondary process */
    unsigned long *status;      /* particle status */
#define PART_ACTIVE       0x0001UL
#define PART_LOST         0x0002UL
    } BEAM;

/* The EMISSION_LOG structure is used for emission logging from cathodes */
typedef struct {
  SDDS_DATASET SDDSout;
  long row, maxRows;
  short active;
} EMISSION_LOG;

/* 
 * Structure CATHODE holds information about generation of new particles from 
 * the cathode
 *
 */

typedef struct {
    double z_position;
    double initial_vz, initial_omega, focal_length;
    double inner_radius, outer_radius;
    double current_density, ref_current_density, temperature, work_function;
    long field_emission;
    double field_emission_beta;
    double start_time, stop_time, time_offset;
    long start_step, stop_step, autophase;
    double number_per_step;
    long step_emission_started;
    double QoMC, Q, electrons_per_macroparticle;
    double inner_r2, outer_r2, stiffness;
    long discretize_radii, correction_interval, zoned_emission;
    long n_bins, n_emitted, spread_over_dt;
    long *distribution, *deficit;
    double *area, *rmin, totalArea;
    long halton_radix_dt, halton_radix_r;
    long halton_ID_dt, halton_ID_r;
    /* for emission profile: */
    long n_prof_pts;
    double *t_prof, *A_prof;
    /* for thermal velocities */
    double sigmaVMaxwellian;
    /* for emission log */
    EMISSION_LOG *emissionLog;
    } CATHODE;


/* The EMITTER structure is used to store information for emission of primary
 * electrons from arbitrary surfaces.
 */

/* Values for emitter orientations */
#define EMITS_RIGHT 1
#define EMITS_UP 2
#define EMITS_LEFT 3
#define EMITS_DOWN 4

typedef struct {
  long materialID;
  short *direction;
  long *iz, *ir, nzr;
  double *number_per_step;
  double initial_v;
  double user_current_density, *current_density, temperature, work_function;
  long field_emission;
  double field_emission_beta;
  double start_time, stop_time, time_offset;
  long start_step, stop_step, autophase;
  long step_emission_started;
  double QoMC, Q, electrons_per_macroparticle;
  long n_emitted, spread_over_dt;
  long halton_radix_dt, halton_radix_r, halton_radix_z;
  long halton_ID_dt, halton_ID_r, halton_ID_z;
  double offsetFactor;
  /* for emission profile: */
  long n_prof_pts;
  double *t_prof, *A_prof;
  /* for thermal velocities */
  double sigmaVMaxwellian;
  /* for emission log */
  EMISSION_LOG *emissionLog;
} EMITTER;

/* The SECONDARY_EMISSION structure is used to store information for
 * emission of secondary electrons.
 */

typedef struct {
  double *energyData, *yieldData, pMin, pMax;
  double emittedMomentum, yieldLimit;
  long nData, verbosity, initialized;
  long maxQueueLength, queueLength;
  short materialID;
  char *inputFile, *logFile;
  SDDS_DATASET SDDSlog;
  /* This is the "queue" of particles that got lost and will generate secondaries */
  double *z, *r, *p, *Q;
  short *generation, *directionCode;
  long *particleID;
} SECONDARY_EMISSION;

/* 
 * The SNAPSHOT structure is used to store information for a beam-snapshot
 * definition.
 * 
 */

typedef struct {
    char *filename;
    double start_time, time_interval;
    long start_step, step_interval;
    long snapshot_index;
    SDDS_TABLE table;
    } SNAPSHOT ;

/* 
 * The SCREEN structure is used to store information about beam-imaging screens
 * that are placed in the beam path
 *
 */

#define N_SCREEN_QUANS 8

typedef struct {
    char *filename;
    double z_position;
    long direction;
    double start_time;
    SDDS_TABLE table;
    double **buffer;    /* (r, pz, pr, pphi, t, q, r0, t0) */
    short *generation;
    int32_t *particleID;
    double *rmin, *rmax, *zmin, *zmax;
    long nb, max_nb;
    } SCREEN;

/*
 * The POINTSP structure is used to store a single element of the AUTOMESH-like
 * input data that defines the metal boundary of a problem.
 *
 */

typedef struct {
    int type, sense;
    double x, y;
    double x0, y0, r, theta, a, b, aOverb;
    double potential;
    short material, ramp_potential;
    } POINTSP;

/*
 * Structure FIELD_SAVING holds information defining desired field saving
 *
 */

typedef struct {
    char *filename;
    double time_interval, start_time;
    double zEzTrigger;
    long step_interval;      /* time_interval/dt */
    long start_step;         /* time_start/dt */
    long save_index;
    long save_before_exiting;
    SDDS_TABLE outTable;
    double EzTriggerLast;
    } FIELD_SAVING;

/*
 * Structure FIELD_OUTPUT holds information defining desired SDDS field output
 *
 */

typedef struct {
    char *filename;
    double time_interval, start_time;
    long step_interval;      /* time_interval/dt */
    long start_step;         /* time_start/dt */
    long output_index, z_interval, r_interval;
    long exclude_imposed_field, separate_imposed_field;
    SDDS_TABLE SDDSout;
    } FIELD_OUTPUT;

/*
 * Structure TEST_BEAM holds information defining the test beam
 * 
 */

typedef struct {
    double start_time;
    long start_step;
    double stop_time;
    long stop_step;
    double initial_z, initial_r;
    double initial_pz, initial_pr, initial_rpphi;
    double z_force;                          /* d/dt(betaz*gamma) */
    double r_force;                          /* d/dt(betar*gamma) */
    double z_start_deceleration;           
    double r_start_deceleration;           
    double electrons_per_macroparticle;    /* 1 macroparticle emitted per time step */
    } TEST_BEAM_DEFN;

/*
 * Function prototypes begin here 
 *
 */

/* prototypes for pic.c */
extern char *compose_filename(char *template, char *root_name);
extern void check_heap(void);

/* prototypes for geometry.c */
void process_geometry_definition(FIELDS *EM_problem, NAMELIST_TEXT *nl_text);

/* prototypes for grid_drawing.c */
void set_up_drawing(int _nx, int _ny, double _xmin, double _xmax, double _ymin, double _ymax, double tolerance_fraction,
    short **_lattice, double **_potential, short **_material);
void starting_point(double x, double y, double potential, short material);
void draw_line(double xi, double yi, double xf, double yf, double potential, short material, double startPotential);
void draw_circle(double xi, double yi, double xc, double yc, double r, double thetaf, double *xf, double *yf, double potential, short material);
void draw_ellipse(double xi, double yi, double xc, double yc, double a, double b, double thetaf, double *xf, double *yf, double potential, short material);
void add_point(double x, double y, double potential, short material);
void add_breakpoint(double x, double y, double potential, short material);
void dump_boundary_points(char *filename, char *rootname);
void dump_dboundary_points(char *filename, char *rootname);
long getBoundaryPoints(long materialID, long **iz, long **ir, short **dir);
void dump_interior_points(char *filename, char *rootname);
void fill_in_drawing(short **lattice, double **psi);
void type_drawing(FILE *fp, char *label, short **lattice, int nx, int ny);
void type_flag_drawing(FILE *fp, char *label, short **_lattice, int nx, int ny);
void trim_drawing(short **lattice);
long out_of_bounds(void);
void set_vacuum_point(double x, double y);
void dumpSurfacePointData(char *rootname, FIELDS *EM);

/* prototypes for trace.c */
void process_trace_request(NAMELIST_TEXT *nltext);
void log_entry(char *routine);
void log_exit(char *routine);

/* flag word for trace mode */
extern long trace_mode;
#define TRACE_ENTRY 1
#define TRACE_HEAP_VERIFY 2
#define TRACEBACK_ON 4

/* prototypes for antenna.c */
void process_antenna_definition(ANTENNA **antenna, long n_antennae, FIELDS *EM_problem, NAMELIST_TEXT *nl_text);
void add_antenna_excitation(FIELDS *EM_problem, ANTENNA *antenna, long n_antennae, double time);

/* prototypes for constant_fields.c */
void set_constant_fields(FIELDS *EM_problem, NAMELIST_TEXT *nl_text);

/* prototypes for field_sampling.c */
void process_sampling_definition(SAMPLE_DEFN **sample_defn, long n_sample_defns, FIELDS *EM_problem, NAMELIST_TEXT *nl_text); 
void sample_fields(SAMPLE_DEFN *sample_defn, long n_sample_defns, FIELDS *EM_problem); 

/* prototypes for integrate.c */
long isqr(long i);
void advance_field_solutions(FIELDS *EM_problem);
void zero_current_arrays(FIELDS *EM_problem);
void multiply_current_arrays(FIELDS *EM_problem, double factor);
void perform_integration(FIELDS *EM_problem, BEAM *beam, CATHODE *cathode, 
                         long cathodes, TEST_BEAM_DEFN *test_beam_defn,
                         ANTENNA *antenna, long n_antennae, RESISTOR *resistor, 
                         SAMPLE_DEFN *sample_defn, long n_sample_defns, 
                         FIELD_SAVING *saving_defn, SNAPSHOT *snapshot_defn, 
                         SCREEN *screen_defn, long n_screen_defns,
                         FIELD_OUTPUT *field_output_defn, SECONDARY_EMISSION *secondaryEmission,
                         long nSecondaryDefns, EMITTER *emitter, long nEmitters, NAMELIST_TEXT *nl_text); 
void smooth_fields(FIELDS *EM_problem, double smoothing_parameter);
void smooth_currents(FIELDS *EM_problem, double smoothing_parameter);
void smooth_currents_SG(FIELDS *EM_problem, long SG_HW);

/* prototypes for field_output.c */
void output_fields(FIELD_OUTPUT *output_defn, FIELDS *EM_problem);
void process_field_output_definition(FIELD_OUTPUT *output_defn, FIELDS *EM_problem, NAMELIST_TEXT *nl_text);

/* prototypes for field_saving.c */
int write_array(char **data, long size_of_element, long n1, long n2, FILE *fp);
void save_fields(FIELD_SAVING *saving_defn, FIELDS *EM_problem, long final_call);
void process_field_saving_definition(FIELD_SAVING *saving_defn, FIELDS *EM_problem, NAMELIST_TEXT *nl_text);

/* prototypes for load_fields.c */
void perform_field_loading(FIELDS *EM_problem, NAMELIST_TEXT *nl_text);
int read_array(char **data, long size_of_element, long n1, long n2, FILE *fp);

/* prototypes for cathode.c */
void process_cathode_definition(FIELDS *EM_problem, CATHODE **cathode, long *cathodes,
                                NAMELIST_TEXT *nl_text);
void emit_electrons_from_cathode(BEAM *beam, CATHODE *cathode, FIELDS *EM_problem, long id);
void set_current_density(CATHODE *cathode, FIELDS *EM_problem);
void addThermalVelocities(BEAM *beam, long ip, double sigmaV);
void logCathodeEmission(EMISSION_LOG *emissionLog, BEAM *beam, long ip, long emitterID);
void setUpEmissionLog(EMISSION_LOG **emissionLog, char *filename);
void finishCathodeEmissionLog(EMISSION_LOG *emissionLog);
double fieldEmissionCurrentDensity(double E, double beta, double wf);

/* prototypes for emitter.c */
void process_emitter_definition(FIELDS *EM_problem, EMITTER **emitter, long *emitters,
                                NAMELIST_TEXT *nl_text);
void emit_electrons_from_emitter(BEAM *beam, EMITTER *emitter, FIELDS *EM_problem);
void set_emitter_current_density(EMITTER *emitter, FIELDS *EM_problem);

/* prototypes for snapshots.c */
void process_snapshot_definition(SNAPSHOT *snapshot_defn, FIELDS *EM_problem, NAMELIST_TEXT *nl_text);
void take_beam_snapshot(SNAPSHOT *snapshot_defn, BEAM *beam, FIELDS *EM_problem);

/* prototypes for advance_particles.c */
void advance_particles(BEAM *beam, FIELDS *EM_problem, TEST_BEAM_DEFN *test_beam_defn, SCREEN *screen_defn, 
                       long n_screens, SDDS_DATASET *SDDSlost, double zForcesStart, 
                       SECONDARY_EMISSION *secondaryEmission, long nSecondaryDefns);
void setParticleVelocities(double *vz, double *vr, double *vphi, double *gamma,
                           double pz, double pr, double r_pphi, double r);
void interpolateFields(double *EzReturn, double *ErReturn, double *EphiReturn,
                       double *BzReturn, double *BrReturn, double *BphiReturn,
                       double z, double r,
                       FIELDS *EM_problem);

/* prototypes for screens.c */
long process_screen_definition(SCREEN **screen_defn, long n_screen_defns, FIELDS *EM_problem,
    NAMELIST_TEXT *nl_text);

/* prototypes for resistor.c */
void process_resistor_definition(RESISTOR *resistor, FIELDS *EM_problem, NAMELIST_TEXT *nl_text);
void add_resistor_currents(FIELDS *EM_problem, RESISTOR *resistor);

/* prototypes for test_beam.c */
void process_testbeam_definition(TEST_BEAM_DEFN *test_beam_defn, FIELDS *EM_problem, 
    NAMELIST_TEXT *nl_text);
void emit_test_particles(BEAM *beam, TEST_BEAM_DEFN *test_beam_defn, FIELDS *EM_problem);

/* prototypes for poisson.c */
long solve_poisson_cyl(double **Psi, double **driving_term, 
    short **flags, long left_bc, long right_bc, 
    long nz, long nr, double dz, double dr, double accuracy, long iteration_limit);

/* prototypes for poisson_correction.c */
void setup_poisson_correction(FIELDS *EM_problem, NAMELIST_TEXT *nl_text);
void add_poisson_correction(FIELDS *EM_problem, BEAM *beam, long space_charge);
void compute_imposed_Efield(FIELDS *EM_problem, long initialize);

/* prototypes for space_charge.c */
void setup_space_charge(FIELDS *EM_problem, NAMELIST_TEXT *nl_text);
void add_particle_currents(FIELDS *EM_problem, BEAM *beam);

/* prototypes for status_pr.c */
void do_status_printouts(FIELDS *EM_problem, BEAM *beam, CATHODE *cathode,
                         long check_divergence, long space_charge,
                         SDDS_DATASET *SDDSout, long *outputRow);
double kinetic_energy(BEAM *beam);
double *field_energy(FIELDS *EM_problem);
void setUpStatusOutputFile(SDDS_DATASET *SDDSout, char *filename, char *rootname);
void setUpLostParticleFile(SDDS_DATASET *SDDSout, char *filename, char *rootname);

/* prototypes for dc_fields.c */
void set_dc_fields(FIELDS *EM_problem, NAMELIST_TEXT *nl_text);

void setup_translation(FIELDS *EM_problem, NAMELIST_TEXT *nl_text);
void translate_problem(BEAM *beam, FIELDS *EM_problem, ANTENNA *antenna, long n_antennae, RESISTOR *resistor);

/* prototypes for load_particles.c */
void perform_particle_loading(FIELDS *EM_problem, BEAM *beam, NAMELIST_TEXT *nl_text);

/* prototypes for solenoid.c and Bsolenoid.c */
void B_ring(double *Brho, double *Bz, double rho, double z, double radius);
void B_solenoid(double *B_rho, double *B_z, double rho, double z, 
               double radius, double zStart, double zEnd, long coils, long symmetry);
void dumpSolenoidFields(char *filename, char *rootname, double **Br, double **Bz, double zmin, double dz, long nz, double dr, 
                        long nr, int32_t count);
void process_solenoid_definition(FIELDS *EM_problem, NAMELIST_TEXT *nl_text);

/* prototypes for onAxisFields.c */
void process_on_axis_field_definition(FIELDS *EM_problem,  NAMELIST_TEXT *nl_text);
void take_derivative(double *dFdz, double *F, double *z, long n_pts, long check_dz, double dzFracLimit);
void addOffAxisExpansionFields(double *Ez, double *Er, double *Bphi, 
                               double *EzMin, double *EzMax, 
                               double z, double r, double t,
                               ON_AXIS_FIELD *onAxisField, long onAxisFields,
                               long profileOnly);

/* prototypes for secondaryEmission.c */
void processSecondaryEmissionDefinition(SECONDARY_EMISSION **secondaryEmission, long *nSecondaryDefns,
                                           NAMELIST_TEXT *nl_text, FIELDS *EM_problem);
void queueSecondaryEmission(SECONDARY_EMISSION *secondaryEmission, long nSecondaryDefns,
                            double p, double z, double r, short directionCode, double Q, short generation,
                            FIELDS *EM_problem, BEAM *beam);
#define GOING_LEFT 1
#define GOING_RIGHT 2
#define GOING_UP 3
#define GOING_DOWN 4
void findLossCoordinates(double *zReturn, double *rReturn, short *directionCode,
                         double z0, double r0, double z1, double r1,
                         FIELDS *EM_problem);
void emitSecondaryElectrons(BEAM *beam, SECONDARY_EMISSION *secondaryEmission,  long nSecondaryDefns,
                            FIELDS *EM_problem);



void setEmitterCurrentDensity(EMITTER *emitter, FIELDS *EM_problem);
