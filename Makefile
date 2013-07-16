# Makefile for SOG code (coupled bio-physical model of deep estuaries)

# Name of the executable we're building
EXEC = SOG

# Compilation
F90 = gfortran
# Don't compile with optimization until the code runs properly without it
# and always revert to -O0 and all checks when adding new features
FFLAGS-DEV = -O0 -g -fimplicit-none -fbounds-check -Wall
# g95 extra flags
#FFLAGS-EXTRA = -ftrace=full
# gfortran extra flags
FFLAGS-EXTRA = -fbacktrace
FFLAGS-PROD = -O3 -fimplicit-none -Wall
LD = $(F90)
LDFLAGS = -o

# Archiver to build objects library for unit tests
AR = ar
ARFLAGS = cr
OBJLIB = libSOG.a
UNITTESTS = tests/unit_tests

# File deleter and flags for cleanup
RM = rm
RMFLAGS = -rf

# Emacs tag generator and flags
ETAGS = etags
ETFLAGS =

# List of objects (order matters)
OBJS = precision_defs.o io_unit_defs.o datetime.o input_processor.o \
fundamental_constants.o malloc.o unit_conversions.o \
grid.o numerics.o forcing.o core_variables.o water_properties.o \
irradiance.o turbulence.o air_sea_fluxes.o carbonate.o fitbottom.o freshwater.o northern_influence.o buoyancy.o \
mixing_layer.o diffusion.o baroclinic_pressure.o upwelling.o \
physics_eqn_builder.o physics_model.o \
NPZD.o biology_eqn_builder.o rungekutta.o biology_model.o \
IMEX_solver.o chemistry_model.o \
profiles_output.o user_output.o timeseries_output.o\
SOG.o

# The executable is the default target that is built by "make"
# It depends on all of the objects which are built from the
# dependencies list by the suffix-based rules below
# The compiler is set to -O0, and lots of checking flags
# (i.e. for development and testing)
$(EXEC): FFLAGS = $(FFLAGS-DEV) $(FFLAGS-EXTRA)
$(EXEC): $(OBJS)
	$(LD) $(OBJS) $(LDFLAGS) $@

# "make tests" creates or updates the library of object files, then
# builds and runs the unit test suite
tests: $(OBJLIB) $(UNITTESTS)/*.f90
	(cd $(UNITTESTS) && make)

$(OBJLIB): FFLAGS = $(FFLAGS-DEV) $(FFLAGS-EXTRA)
$(OBJLIB): $(OBJS)
	$(AR) $(ARFLAGS) $(OBJLIB) $(OBJS)

# "make tags" builds a list of tags for subroutines, modules, etc. to
# make source file navigation in Emacs easier.  By convention the tags
# are stored in a file called TAGS
.PHONY: tags
tags:
	$(ETAGS) $(ETFLAGS) *.f90

# "make clean" deletes all object, module files, any core dumps, the
# executable, and the object library
.PHONY: clean
clean:
	$(RM) $(RMFLAGS) *.o *.mod core $(EXEC) $(OBJLIB)

# "make changelog" builds or updates a GNU-style ChangeLog file
# from the CVS log messages
.PHONY: changelog
changelog:
	$(CVS2CL) $(CLFLAGS)

# "make SOG-dev" does a clean build with the compiler flags set to
# -O0, and lots of checking (i.e. appropriate for development and testing)
$(EXEC)-dev: FFLAGS = $(FFLAGS-DEV) $(FFLAGS-EXTRA)
$(EXEC)-dev: clean $(OBJS)
	$(LD) $(OBJS) $(LDFLAGS) $(EXEC)

# "make SOG-prod" does a clean build with the compiler flags set
# for faster execution.  *** Don't use this until you've sure the code
# is working in development mode, and you've compared the results of a
# few production and development builds to ensure that the code is stable
# with optimization enabled. ***  Consider yourself warned!
$(EXEC)-prod: FFLAGS = $(FFLAGS-PROD)
$(EXEC)-prod: clean $(OBJS)
	$(LD) $(OBJS) $(LDFLAGS) $(EXEC)

# Suffix-based compilation rules
.SUFFIXES: 		# Delete built-in implicit suffix-based rules
.SUFFIXES: .f90 .o

.f90.o:
	$(F90) -c $(FFLAGS) $<
