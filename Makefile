# $Id$
# $Source$

# Makefile for SoG code (bio-physical model of the Strait of Georgia)
# Configured for g95 on Mac OS/X and for pgf90 on coho at UBC EOS

# Name of the executable we're building
EXEC = SOG

# The Tool List
# Fortran 90 compiler, linker/loader, and their flags
ifeq "$(HOST)" "coho"
  # Use pgf90 compiler on coho
  F90 = pgf90
  # Don't compile with optimization until the code runs properly without it
  # and always revert to -O0 and all checks when adding new features 
  FFLAGS = -O0 -g -Mbounds -Mdclchk -Mstandard -Minform=warn
  LD = pgf90
  LDFLAGS = -o
else
  # Use g95 compile (intended for iBooks, but might work on other platforms)
  F90 = g95
  # Don't compile with optimization until the code runs properly without it
  # and always revert to -O0 and all checks when adding new features
  FFLAGS = -O0 -g -fimplicit-none -fbounds-check -ftrace=full -Wall
  LD = g95
  LDFLAGS = -o
endif
# File deleter and flags for cleanup
RM = rm
RMFLAGS = -rf
# Emacs tag generator and flags
ETAGS = etags
ETFLAGS = 
# ChangeLog file generator
CVS2CL = ../cvs2cl/cvs2cl.pl
CLFLAGS = -rt

# List of objects (order matters)
# timeseries_output.o
OBJS = precision_defs.o io_unit_defs.o malloc.o unit_conversions.o \
datetime.o grid.o \
\
mean_param.o declarations.o surface_forcing.o \
\
water_properties.o find_upwell.o diffusion.o fitbottom.o \
getpar.o \
\
define_flux.o \
initial_sog.o IMEX_constants.o write_open.o           \
allocate1.o read_sog.o find_wind.o 	\
stability.o \
\
allocate2.o Julian_day.o initialize.o	\
smoothdata.o coefficients.o define_grid.o define_sog.o	\
alpha_sub.o density_sub.o irradiance_sog.o 	\
buoyancy.o surface_flux_sog.o fun_constants.o		\
ND_flux_profile.o vel_scales.o convection_scales.o		\
shear_diff.o double_diff.o interior_match.o interior_match2.o		\
shape_parameters.o modify_K.o def_gamma.o 	\
Coriolis_and_pg.o 		\
matrix_A.o scalar_H.o U_H.o def_v_t_sog.o TRIDAG.o	\
define_Ri_b_sog.o ML_height_sog.o define_hm_sog.o pdf.o 		\
write_physical_sog.o new_year.o SOG.o define_PZ.o		\
odeint.o reaction_p_sog.o P_H.o find_new.o			\
fit.o polint.o gammq.o gser.o gcf.o	\
gammln.o allocate3.o reaction.o derivs_noflag.o derivs_sog.o rkqs.o     \
rkck.o find_PON.o advection.o

# The executable is the default target that is built by "make"
# It depends on all of the objects which are built from the
# dependencies list by the suffix-based rules below
$(EXEC): $(OBJS)
	$(LD) $(OBJS) $(LDFLAGS) $@

# "make tags" builds a list of tags for subroutines, modules, etc. to
# make source file navigation in Emacs easier.  By convention the tags
# are stored in a file called TAGS
.PHONY: tags
tags:
	$(ETAGS) $(ETFLAGS) *.f90 *.f

# "make clean" deletes all object, and module files,
# any core dump file, and the executable
.PHONY: clean
clean:
	$(RM) $(RMFLAGS) *.o *.mod core $(EXEC)

# "make changelog" builds or updates a GNU-style ChangeLog file
# from the CVS log messages
.PHONY: changelog
changelog:
	$(CVS2CL) $(CLFLAGS)

# Suffix-based compilation rules
.SUFFIXES: 		# Delete built-in implicit suffix-based rules
.SUFFIXES: .f90 .f .o

.f90.o:
	$(F90) -c $(FFLAGS) $<
.f.o: 
	$(F90) -c $(FFLAGS) $<

# end of file
