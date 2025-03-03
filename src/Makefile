SDDS_REPO = $(wildcard ../../SDDS)
ifeq ($(SDDS_REPO),)
  $(error SDDS source code not found. Run 'git clone https://github.com/rtsoliday/SDDS.git' next to the spiffe repository)
endif

include ../Makefile.rules

PROD = spiffe

spiffe_SRC = advance_particles.c antenna.c cathode.c constant_fields.c \
		field_sampling.c field_saving.c  field_output.c geometry.c \
		grid_drawing.c integrate.c load_fields.c load_particles.c \
		poisson.c poisson_correction.c resistor.c screens.c snapshots.c \
		test_beam.c translate.c space_charge.c spiffe.c status_pr.c \
		Bsolenoid.c solenoid.c onAxisFields.c secondaryEmission.c emitter.c

SVN_VERSION=GIT_$(shell git log --oneline | wc -l | sed 's/^[[:space:]]*//')

CFLAGS += -I$(SDDS_REPO)/include -I$(OBJ_DIR) -DSVN_VERSION=\"$(SVN_VERSION)\"

ifeq ($(OS), Linux)
  CFLAGS += -Wno-unused-variable
  LDFLAGS := -L$(SDDS_REPO)/lib/$(OS)-$(ARCH) -fopenmp $(LDFLAGS)
  PROD_SYS_LIBS := $(LZMA_LIB) $(GSL_LIB) $(GSLCBLAS_LIB) $(Z_LIB) $(PROD_SYS_LIBS)
  PROD_LIBS_SDDS = -lmdbcommon -lmatlib -lSDDS1 -lnamelist -lrpnlib -lmdbmth -lmdblib
endif

ifeq ($(OS), Darwin)
  CFLAGS += -Wno-unused-variable
  LDFLAGS := -L$(SDDS_REPO)/lib/$(OS)-$(ARCH) $(LDFLAGS)
  PROD_SYS_LIBS := $(LZMA_LIB) $(GSL_LIB) $(GSLCBLAS_LIB) $(Z_LIB) $(PROD_SYS_LIBS)
  PROD_LIBS_SDDS = -lmdbcommon -lmatlib -lSDDS1 -lnamelist -lrpnlib -lmdbmth -lmdblib
endif

ifeq ($(OS), Windows)
  CFLAGS += -I$(SDDS_REPO)/lzma -wd4244
  PROD_LIBS += gsl.lib gslcblas.lib
  PROD_LIBS_SDDS = mdbcommon.lib matlib.lib SDDS1.lib namelist.lib rpnlib.lib mdbmth.lib mdblib.lib lzma.lib z.lib
  LIB_LINK_DIRS += -LIBPATH:$(SDDS_REPO)/lib/$(OS)-$(ARCH)
endif

include ../Makefile.build

$(OBJ_DIR)/spiffe$(EXEEXT): $(spiffe_OBJS) $(PROD_DEPS)
	$(LINKEXE) $(OUTPUTEXE) $(spiffe_OBJS) $(LDFLAGS) $(LIB_LINK_DIRS) $(PROD_LIBS) $(PROD_LIBS_SDDS) $(PROD_SYS_LIBS)
	cp -f $@ $(BIN_DIR)/

$(OBJ_DIR)/antenna.$(OBJEXT): antenna.c $(OBJ_DIR)/antenna.h

$(OBJ_DIR)/antenna.h: antenna.nl
	../../SDDS/bin/$(OS)-$(ARCH)/nlpp -suppressSummaryVariables antenna.nl $(OBJ_DIR)/antenna.h

$(OBJ_DIR)/cathode.$(OBJEXT): cathode.c $(OBJ_DIR)/cathode.h

$(OBJ_DIR)/cathode.h: cathode.nl
	../../SDDS/bin/$(OS)-$(ARCH)/nlpp -suppressSummaryVariables cathode.nl $(OBJ_DIR)/cathode.h

$(OBJ_DIR)/constant_fields.$(OBJEXT): constant_fields.c $(OBJ_DIR)/constant_fields.h

$(OBJ_DIR)/constant_fields.h: constant_fields.nl
	../../SDDS/bin/$(OS)-$(ARCH)/nlpp -suppressSummaryVariables constant_fields.nl $(OBJ_DIR)/constant_fields.h

$(OBJ_DIR)/field_output.$(OBJEXT): field_output.c $(OBJ_DIR)/field_output.h

$(OBJ_DIR)/field_output.h: field_output.nl
	../../SDDS/bin/$(OS)-$(ARCH)/nlpp -suppressSummaryVariables field_output.nl $(OBJ_DIR)/field_output.h

$(OBJ_DIR)/field_sampling.$(OBJEXT): field_sampling.c $(OBJ_DIR)/field_sampling.h

$(OBJ_DIR)/field_sampling.h: field_sampling.nl
	../../SDDS/bin/$(OS)-$(ARCH)/nlpp -suppressSummaryVariables field_sampling.nl $(OBJ_DIR)/field_sampling.h

$(OBJ_DIR)/field_saving.$(OBJEXT): field_saving.c $(OBJ_DIR)/field_saving.h

$(OBJ_DIR)/field_saving.h: field_saving.nl
	../../SDDS/bin/$(OS)-$(ARCH)/nlpp -suppressSummaryVariables field_saving.nl $(OBJ_DIR)/field_saving.h

$(OBJ_DIR)/geometry.$(OBJEXT): geometry.c $(OBJ_DIR)/geometry.h

$(OBJ_DIR)/geometry.h: geometry.nl
	../../SDDS/bin/$(OS)-$(ARCH)/nlpp geometry.nl $(OBJ_DIR)/geometry.h

$(OBJ_DIR)/integrate.$(OBJEXT): integrate.c $(OBJ_DIR)/integrate.h

$(OBJ_DIR)/integrate.h: integrate.nl
	../../SDDS/bin/$(OS)-$(ARCH)/nlpp -suppressSummaryVariables integrate.nl $(OBJ_DIR)/integrate.h

$(OBJ_DIR)/load_fields.$(OBJEXT): load_fields.c $(OBJ_DIR)/load_fields.h

$(OBJ_DIR)/load_fields.h: load_fields.nl
	../../SDDS/bin/$(OS)-$(ARCH)/nlpp -suppressSummaryVariables load_fields.nl $(OBJ_DIR)/load_fields.h

$(OBJ_DIR)/load_particles.$(OBJEXT): load_particles.c $(OBJ_DIR)/load_particles.h

$(OBJ_DIR)/load_particles.h: load_particles.nl
	../../SDDS/bin/$(OS)-$(ARCH)/nlpp -suppressSummaryVariables load_particles.nl $(OBJ_DIR)/load_particles.h

$(OBJ_DIR)/onAxisFields.$(OBJEXT): onAxisFields.c $(OBJ_DIR)/onAxisFields.h

$(OBJ_DIR)/onAxisFields.h: onAxisFields.nl
	../../SDDS/bin/$(OS)-$(ARCH)/nlpp -suppressSummaryVariables onAxisFields.nl $(OBJ_DIR)/onAxisFields.h

$(OBJ_DIR)/poisson_correction.$(OBJEXT): poisson_correction.c $(OBJ_DIR)/poisson_correction.h

$(OBJ_DIR)/poisson_correction.h: poisson_correction.nl
	../../SDDS/bin/$(OS)-$(ARCH)/nlpp -suppressSummaryVariables poisson_correction.nl $(OBJ_DIR)/poisson_correction.h

$(OBJ_DIR)/resistor.$(OBJEXT): resistor.c $(OBJ_DIR)/resistor.h

$(OBJ_DIR)/resistor.h: resistor.nl
	../../SDDS/bin/$(OS)-$(ARCH)/nlpp -suppressSummaryVariables resistor.nl $(OBJ_DIR)/resistor.h

$(OBJ_DIR)/screens.$(OBJEXT): screens.c $(OBJ_DIR)/screens.h

$(OBJ_DIR)/screens.h: screens.nl
	../../SDDS/bin/$(OS)-$(ARCH)/nlpp -suppressSummaryVariables screens.nl $(OBJ_DIR)/screens.h

$(OBJ_DIR)/snapshots.$(OBJEXT): snapshots.c $(OBJ_DIR)/snapshots.h

$(OBJ_DIR)/snapshots.h: snapshots.nl
	../../SDDS/bin/$(OS)-$(ARCH)/nlpp -suppressSummaryVariables snapshots.nl $(OBJ_DIR)/snapshots.h

$(OBJ_DIR)/space_charge.$(OBJEXT): space_charge.c $(OBJ_DIR)/space_charge.h

$(OBJ_DIR)/space_charge.h: space_charge.nl
	../../SDDS/bin/$(OS)-$(ARCH)/nlpp -suppressSummaryVariables space_charge.nl $(OBJ_DIR)/space_charge.h

$(OBJ_DIR)/test_beam.$(OBJEXT): test_beam.c $(OBJ_DIR)/test_beam.h

$(OBJ_DIR)/test_beam.h: test_beam.nl
	../../SDDS/bin/$(OS)-$(ARCH)/nlpp -suppressSummaryVariables test_beam.nl $(OBJ_DIR)/test_beam.h

$(OBJ_DIR)/translate.$(OBJEXT): translate.c $(OBJ_DIR)/translate.h

$(OBJ_DIR)/translate.h: translate.nl
	../../SDDS/bin/$(OS)-$(ARCH)/nlpp -suppressSummaryVariables translate.nl $(OBJ_DIR)/translate.h

$(OBJ_DIR)/solenoid.$(OBJEXT): solenoid.c $(OBJ_DIR)/solenoid.h

$(OBJ_DIR)/solenoid.h: solenoid.nl
	../../SDDS/bin/$(OS)-$(ARCH)/nlpp -suppressSummaryVariables solenoid.nl $(OBJ_DIR)/solenoid.h

$(OBJ_DIR)/secondaryEmission.$(OBJEXT): secondaryEmission.c $(OBJ_DIR)/secondaryEmission.h

$(OBJ_DIR)/secondaryEmission.h: secondaryEmission.nl
	../../SDDS/bin/$(OS)-$(ARCH)/nlpp -suppressSummaryVariables secondaryEmission.nl $(OBJ_DIR)/secondaryEmission.h

$(OBJ_DIR)/emitter.$(OBJEXT): emitter.c $(OBJ_DIR)/emitter.h

$(OBJ_DIR)/emitter.h: emitter.nl
	../../SDDS/bin/$(OS)-$(ARCH)/nlpp -suppressSummaryVariables emitter.nl $(OBJ_DIR)/emitter.h

