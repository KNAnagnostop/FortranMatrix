#====================================================================================================================
# $@ target, $< first dependency, $^ all dependencies, 
# $* "stem" of target filename e.g. %.f: %.F will five foo for foo.f
# @command executes command without echoing it; make -s (silent make); make --no-print-directory (no echoing of dirs)
#====================================================================================================================
.PHONY: all clean mods new debug test testvec
#======================================================================
SYSTEM:=$(shell uname)
SHELL  =  /bin/bash
MATRIX = $(HOME)/lib/Matrix/
#======================================================================
OBJS   = 
MODS   =  $(subst .o,.mod,$(OBJS))
FC     =  gfortran
IFLAGS = -I$(MATRIX)/lib/gnu/
FFLAGS = -Ofast -march=native -mtune=native -ffast-math -funroll-loops
DFLAGS = -g -ffpe-trap=zero,invalid,overflow,underflow,denormal  -fdump-core -fbacktrace
LIBS   = -L$(MATRIX)/lib/gnu/ -lMatrix -llapack
#----------------------------------------------------------------------
IFC     =  ifort
IFFLAGS =  -O3 -heap-arrays -xHost  -parallel      
ILIBS   =  -mkl=parallel
INTEL   =  /opt/intel/oneapi/setvars.sh
#======================================================================

all:  testpauli


test:    
	$(FC)   $(IFLAGS) $(DFLAGS) test_basic.f90       $(LIBS) -o t

testvec:    
	$(FC)   $(IFLAGS) $(DFLAGS) test_vec.f90         $(LIBS) -o t

testasym:    
	$(FC)   $(IFLAGS) $(DFLAGS) test_asym.f90        $(LIBS) -o t

testvectorize:    
	$(FC)   $(IFLAGS) $(DFLAGS) test_vectorize.f90   $(LIBS) -o t

testpauli:
	$(FC)   $(IFLAGS) $(DFLAGS) test_pauli.f90       $(LIBS) -o t

new:
	@make -s                   clean
	@make --no-print-directory

debug:
	@make -s                   clean
	@make --no-print-directory FFLAGS="$(DFLAGS)"

intel:
	. $(INTEL); make --no-print-directory FC=$(IFC) FFLAGS="$(IFFLAGS) $(ILIBS)"
mods:
	@make --no-print-directory mods $(MODS)

%.mod:  %.f90
	$(FC) -c $(FFLAGS) $<

%.o:    %.f90
	$(FC) -c $(FFLAGS) $<  -o $@

clean:
	$(RM) *.o *.mod t
