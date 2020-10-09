F90          = gfortran
F90LINKER    = gfortran

#F90 = ifort
#F90LINKER = ifort

# for gfortran Compiler
#======================
#DEFS      =
#FFLAGS   = -fdefault-real-8 -fdefault-double-8 -fbounds-check -fbacktrace -g
FFLAGS   	= -Og -pipe -g -fbacktrace -Wall -Wextra -pedantic -fcheck=all -Wconversion -fbounds-check
#FFLAGS   = -O3 -pipe
#INCLUDES  =
#LFLAGS    = $(FFLAGS)

# for ifort Compiler
#===================
DEFS      =
#FFLAGS   = -O0 -check all -C -warn -qopenmp -g -traceback -fpp -prec-div -fp-model source -fpe0
#FFLAGS    = -O3 -ipo -xHost -qopenmp -fpp -fp-model source
INCLUDES  =
LFLAGS    = $(FFLAGS)

OBJECTS = \
  k_Rosseland_mod.o \
	IC_mod.f90

# executable statement
EXECS  = IC

.SUFFIXES : .o .f .F .f90

default: IC

IC:  $(OBJECTS)
	 $(F90LINKER) $(LFLAGS) $(OBJECTS) -o $(EXECS)

clean:
	rm -f *.o

realclean:
	rm -f *.o *~ *.mod *__genmod.f90 $(EXECS)

.f.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<
.F.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<
.f90.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<
