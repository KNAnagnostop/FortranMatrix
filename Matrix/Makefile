#====================================================================================================================
# $@ target, $< first dependency, $^ all dependencies, 
# $* "stem" of target filename e.g. %.f: %.F will five foo for foo.f
# @command executes command without echoing it; make -s (silent make); make --no-print-directory (no echoing of dirs)
#====================================================================================================================
# export OMP_NUM_THREADS=8; export MKL_NUM_THREADS=8
#======================================================================
.PHONY: all clean mods new debug d n lib intel funclist garis pack
#======================================================================
SYSTEM :=$(shell uname)
SHELL   =  /bin/bash
#======================================================================
OBJS    =  matrix_mod_common.o matrix_mod_array.o matrix_mod_matrix.o tensorprod_mod.o timer_mod.o array_mod.o
MODS    =  $(subst .o,.mod,$(OBJS))
SRCS    =  $(subst .o,.f90,$(OBJS)) matrix_mod_matrix_matrix.f90 matrix_mod_matrix_matrixClass.f90 matrix_mod_matrix_dmatrix.f90 matrix_mod_matrix_vector.f90 matrix_mod_matrix_dvector.f90  matrix_mod_array_lapack.f90 tensorprod_mod.f90 timer_mod.f90
#----------------------------------------------------------------------
FC      =  gfortran  # caf,ifort,gfortran
FFLAGS  =  -Ofast -march=native -mtune=native -ffast-math -funroll-loops
LIBS    =  -llapack -lblas 
#----------------------------------------------------------------------
IFC     =  ifort
IFFLAGS =  -O3 -heap-arrays -xHost  -parallel      
INTELV  =  . /opt/intel/oneapi/setvars.sh #/opt/intel/bin/compilervars.sh intel64;
IMKL    =  -mkl=parallel                  # version 2022:  -qmkl=parallel
ILIBS   =  $(IMKL)
#======================================================================
# modules to load on ARIS:
ARMODULES =   intelmpi/2018 intel/18     # check if you can load intel/19, there is not intelmpi/2019 yet
# compile with gfortran, load mkl
GARMODULES=   gnu/8 intel/18 intelmpi/2018 
GAFC      =   gfortran
GAFFLAGS  =   -Ofast -march=native -mtune=native -ffast-math -funroll-loops -fdefault-integer-8  -m64 -I"${MKLROOT}/include" 
GALIBS    =   -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_ilp64 -lmkl_tbb_thread -lmkl_core -ltbb -lstdc++ -lpthread -lm -ldl
# Compile with PGI fortran:
PGMODULES = pgi/19.4
PGFC      = pgfortran
PGFFLAGS  = -O4   -Mvect=altcode,fuse,gather,simd,256,tile -m64 -Mconcur=allcores,assoc -Mpropcond -Minfo=all # -tp=haswell # -ta=nvidia -acclibs -acc=autopar # for OpenACC
PGLIBS    = 
#======================================================================

all: 
	make gnu
	make intel

#$(OBJS): $(MODS)

test:   $(OBJS) test.o  
	$(FC)   $(FFLAGS) $^  -o $@	$(LIBS)

testa:   $(OBJS) testa.o  
	$(FC)   $(FFLAGS) $^  -o $@	$(LIBS)

array_benchmark:  $(OBJS) array_benchmark.o
	$(FC)   $(FFLAGS) $^  -o $@	$(LIBS)

matrix_benchmark: $(OBJS) matrix_benchmark.o  
	$(FC)   $(FFLAGS) $^  -o $@	$(LIBS)

matrix_benchmark_mymatrix: $(OBJS) matrix_benchmark_mymatrix.o  
	$(FC)   $(FFLAGS) $^  -o $@	$(LIBS)

# Include files:
matrix_mod_matrix.o: matrix_mod_matrix.f90 matrix_mod_matrix_matrix.f90 matrix_mod_matrix_matrixClass.f90 matrix_mod_matrix_dmatrix.f90 matrix_mod_matrix_vector.f90 matrix_mod_matrix_dvector.f90
matrix_mod_array.o : matrix_mod_array.f90  matrix_mod_array_lapack.f90  
# all target, using intel compiler
i:
	set +u;$(INTELV); make --no-print-directory FC="$(IFC)" FFLAGS="$(IFFLAGS)" LIBS="$(ILIBS)"

intel:
	make clean
	set +u;$(INTELV); make --no-print-directory FC="$(IFC)" FFLAGS="$(IFFLAGS)" LIBS="$(ILIBS)"  mods
	ar rcv libMatrix.a $(OBJS)
	cp libMatrix.a *.mod lib/intel

string:
	make IMKL="-qmkl=parallel" intel

aris:
	make clean
	module purge;module load $(ARMODULES);           make --no-print-directory FC="$(IFC)" FFLAGS="$(IFFLAGS)" LIBS="$(ILIBS)"  mods
	ar rcv libMatrix.a $(OBJS)
	cp libMatrix.a *.mod lib/aris

#OBJS    =  matrix_mod_common.o matrix_mod_array.o matrix_mod_matrix.o array_mod.o
pgiaris:
	make clean
	module purge;module load $(PGMODULES);           make --no-print-directory FC="$(PGFC)" FFLAGS="$(PGFFLAGS)" LIBS="$(PGLIBS)"  mods;
	ar rcv libMatrix.a $(OBJS)
	cp libMatrix.a *.mod lib/pgiaris


garis:
	make clean
	module purge;module load $(GARMODULES);          make --no-print-directory FC="$(GAFC)" FFLAGS="$(GAFFLAGS)" LIBS="$(GALIBS)"  mods
	ar rcv libMatrix.a $(OBJS)
	cp libMatrix.a *.mod lib/garis

gnu:
	make clean
	make --no-print-directory FC="$(FC)" FFLAGS="$(FFLAGS)" LIBS="$(LIBS)"  mods
	ar rcv libMatrix.a $(OBJS)
	cp libMatrix.a *.mod lib/gnu

n:      new
new:
	@make -s                   clean
	@make --no-print-directory


d:      debug
debug:
	@make -s                   clean
	@make --no-print-directory FFLAGS=-g

mods:
	@make --no-print-directory $(OBJS)

%.mod: %.f90
	$(FC) -c $(FFLAGS) $<

%.o:   %.f90
	$(FC) -c $(FFLAGS) $<  -o $@

# -I/path/to/_.mod/files/  -L/path/to/libMatrix/                                                   (finclude/ is the default name of directories with *.mod files)
# Use of library: Compile your program with: gfortran myprog.f90 -I./Matrix -L./Matrix -lMatrix
lib:
	@make  --no-print-directory mods
	ar rcv libMatrix.a $(OBJS)

funclist:
	@awk '/:\.:\.:/{print substr($$0,0,138)}' $(SRCS) >  matrix_proclist;a2ps -1 --chars-per-line=138 -r -o - matrix_proclist | ps2pdf - > matrix_proclist.pdf

clean:
	$(RM) test test_benchmark test.dat *.o *.mod core libMatrix.a fort.[0-9][0-9] Matrix.tgz

pack:
	@./pack
