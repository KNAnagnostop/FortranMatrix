#====================================================================================================================
# $@ target, $< first dependency, $^ all dependencies, 
# $* "stem" of target filename e.g. %.f: %.F will five foo for foo.f
# @command executes command without echoing it; make -s (silent make); make --no-print-directory (no echoing of dirs)
#====================================================================================================================
# export OMP_NUM_THREADS=8; export MKL_NUM_THREADS=8
#======================================================================
.PHONY: all clean mods new debug d n lib intel funclist aris garis pack pfapack pfapackintel pfapackaris  pfapackgaris libs aris nvaris
#======================================================================
SYSTEM :=$(shell uname)
SHELL   =  /bin/bash
#======================================================================
OBJS    =  matrix_mod_common.o matrix_mod_array.o matrix_mod_matrix.o tensorprod_mod.o timer_mod.o array_mod.o
MODS    =  $(subst .o,.mod,$(OBJS))
SRCS    =  $(subst .o,.f90,$(OBJS)) matrix_mod_matrix_matrix.f90 matrix_mod_matrix_matrixClass.f90 matrix_mod_matrix_dmatrix.f90 matrix_mod_matrix_vector.f90 matrix_mod_matrix_dvector.f90  matrix_mod_array_lapack.f90 tensorprod_mod.f90 timer_mod.f90
#----------------------------------------------------------------------
FC      =  gfortran  # caf,ifort,gfortran
FFLAGS  = -O3 -mtune=native -funroll-loops -fopenmp -pthread             #  (bug in sort routines):-march=native maybe dangerous: -floop-parallelize-all 
#FFLAGS =  -Ofast -mtune=native -funroll-loops -fopenmp -pthread  -fallow-argument-mismatch #   #  (bug in sort routines):-march=native maybe dangerous: -floop-parallelize-all 
LIBS    =  -llapack -lblas 
#----------------------------------------------------------------------
IFC     =  ifx
IFFLAGS =  -O3 -xHost -qopenmp -qmkl=parallel   
INTELV  =  . /opt/intel/oneapi/setvars.sh #/opt/intel/bin/compilervars.sh intel64;
IMKL    =  -qmkl=parallel                 # version 2022:  -qmkl=parallel
AMKL    =   -mkl=parallel                 # ARIS has older intel
ILIBS   =  $(IMKL)
ALIBS   =  $(AMKL)
#======================================================================
# modules to load on ARIS 2025:
# module load aocc/5.0.0 openmpi/4.1.8/aocc aocl/5.1.0; mpifort -showme:command  :: gives the flang compiler.
#----------------------------------------------------------------------
# compile: module purge; module load aocc/5.0.0 openmpi/4.1.8/aocc aocl/5.1.0;make aris
#----------------------------------------------------------------------
ARMODULES   =   aocc/5.0.0 openmpi/4.1.8/aocc aocl/5.1.0
ARFC        =   flang   # flang is the underlying compiler
ARFFLAGS    =   -O3 -march=znver3 -mtune=znver3 -fopenmp -funroll-loops    -flto -fuse-ld=lld # -fuse-ld=lld for correct linking
ARLIBS      =   -lscalapack -lflame -lblis-mt -laoclutils -lstdc++
ARAR        = llvm-ar
ARRANLIB    = llvm-ranlib
ARLD        = ld.lld
#----------------------------------------------------------------------
# Aris MKL using gfortran:
MKARMODULES =  gnu/13 oneAPI/2025.3.0 mpi/2021.17 mkl/2025.3
MKARFC      =  gfortran    # nvfortran is  the underlying compiler
MKARFFLAGS  =  -O3 -march=znver3 -mtune=znver3 -funroll-loops -fopenmp -fallow-argument-mismatch
MKARLIBS    =  -L${MKLROOT}/lib -lmkl_scalapack_lp64 -Wl,-no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -lgomp -lpthread -lm -ldl
#----------------------------------------------------------------------
# Nvidia compiler:
#----------------------------------------------------------------------
# compile: module purge; module load NVHPC_SDK/26.5 nvhpc-hpcx-cuda13/26.5; make nvaris
#----------------------------------------------------------------------
NVARMODULES =  NVHPC_SDK/26.5 nvhpc-hpcx-cuda13/26.5
NVARFC      =  nvfortran    # nvfortran is  the underlying compiler
NVARFFLAGS  =  -O3 -tp=zen3 -Munroll -Mnostack_arrays -mp
NVARLIBS    =   -Mscalapack -llapack -lblas
#======================================================================

#----------------------------------------------------------------------
# Nvidia compiler on ICARUS:
#----------------------------------------------------------------------
ICMODULES =  nvhpc-hpcx-cuda13-openmpi5/26.9
ICFC      =  nvfortran    # nvfortran is  the underlying compiler
#ICFFLAGS =  -fast -O3 -Mvect -Minline -gpu=mem:unified -stdpar=gpu -mp=gpu -Mnvpl -Mscalapack -cudalib=cublas -nvmalloc -Minfo
ICFFLAGS  =  -fast -O4  -Mnostack_arrays -Minline -gpu=mem:unified -stdpar=gpu -mp=gpu # -Minfo
ICLIBS    =  
#======================================================================

#----------------------------------------------------------------------
# Nvidia compiler on ICARUS: no GPU
#----------------------------------------------------------------------
ICPUMODULES =  nvhpc-hpcx-cuda13-openmpi5/26.5
ICPUFC      =  nvfortran    # nvfortran is  the underlying compiler
ICPUFFLAGS  =  -fast -O3 -Mvect -Minline -stdpar=multicore -Mconcur -mp -Mnvpl -Mscalapack -nvmalloc -Minfo
ICPULIBS    =  
#======================================================================

all: 
	make gnu
	make intel

#$(OBJS): $(MODS)
libs:
	$(AR) rc libMatrix.a $(OBJS)

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

pfapack:
	cd pfapack/fortran;cp makefile.gnu makefile;make clean;make; cp *.mod lib*.a ../../lib/gnu/; cp *.mod ../../; $(AR) rc ../../libMatrix.a *.o;$(AR) rc ../../lib/gnu/libMatrix.a *.o
#	cd Pfapack/fortran;cp make.gnu make.inc;make clean;make; cp *.mod lib*.a ../../lib/gnu/; ar rcv ../../lib/gnu/libMatrix.a *.o

# Include files:
matrix_mod_matrix.o: matrix_mod_matrix.f90 matrix_mod_matrix_matrix.f90 matrix_mod_matrix_matrixClass.f90 matrix_mod_matrix_dmatrix.f90 matrix_mod_matrix_vector.f90 matrix_mod_matrix_dvector.f90
matrix_mod_array.o : matrix_mod_array.f90  matrix_mod_array_lapack.f90  

#-----------------------------------------------------------------------------------------------------------------------------
# all target, using intel compiler
i:      
	make pfapackintel
	set +u;type ifx &> /dev/null || $(INTELV); make --no-print-directory FC="$(IFC)" FFLAGS="$(IFFLAGS)" LIBS="$(ILIBS)"

intel:  
	make clean
	make pfapackintel
	set +u;type ifx &> /dev/null || $(INTELV); make --no-print-directory FC="$(IFC)" FFLAGS="$(IFFLAGS)" LIBS="$(ILIBS)"  mods
	$(AR) rc libMatrix.a $(OBJS)
	cp libMatrix.a *.mod lib/intel

pfapackintel:
	cd pfapack/fortran;set +u;type ifx &> /dev/null || $(INTELV);cp makefile.intel makefile;make clean;make; cp *.mod lib*.a ../../;cp *.mod lib*.a ../../lib/intel/; $(AR) rc ../../lib/intel/libMatrix.a *.o; $(AR) rc ../../libMatrix.a *.o


#-----------------------------------------------------------------------------------------------------------------------------
aris:   
	make clean
	module purge;module load $(ARMODULES); make  AR="$(ARAR)" RANLIB="$(ARRANLIB)" LD="$(ARLD)" pfapackaris
	module purge;module load $(ARMODULES); make --no-print-directory FC="$(ARFC)" FFLAGS="$(ARFFLAGS)" LIBS="$(ARLIBS)" AR="$(ARAR)" RANLIB="$(ARRANLIB)" LD="$(ARLD)"  mods
	module purge;module load $(ARMODULES); make  AR="$(ARAR)" RANLIB="$(ARRANLIB)" LD="$(ARLD)" libs
	cp libMatrix.a *.mod lib/aris

pfapackaris:
	cd pfapack/fortran;module purge;module load $(ARMODULES);cp makefile.aris makefile;make clean;make; cp *.mod lib*.a ../../lib/aris/; cp *.mod ../../; $(ARAR) rc ../../lib/aris/libMatrix.a *.o; $(ARAR) rc ../../libMatrix.a *.o

#-----------------------------------------------------------------------------------------------------------------------------
# MKL+GNU on ARIS:
arismkl:   
	make clean
	make pfapackarismkl
	module purge;module load $(MKARMODULES); make --no-print-directory FC="$(MKARFC)" FFLAGS="$(MKARFFLAGS)" LIBS="$(MKARLIBS)"  mods
	make  libs
	cp libMatrix.a *.mod lib/arismkl

pfapackarismkl:
	cd pfapack/fortran;module purge;module load $(MKARMODULES);cp makefile.arismkl makefile;make clean;make; cp *.mod lib*.a ../../lib/arismkl/; cp *.mod ../../; $(AR) rc ../../lib/arismkl/libMatrix.a *.o; $(AR) rc ../../libMatrix.a *.o

#-----------------------------------------------------------------------------------------------------------------------------
# NVIDIA on ARIS:
nvaris:   
	make clean
	make pfapacknvaris
	module purge;module load $(NVARMODULES); make --no-print-directory FC="$(NVARFC)" FFLAGS="$(NVARFFLAGS)" LIBS="$(NVARLIBS)"  mods
	make  libs
	cp libMatrix.a *.mod lib/nvaris

pfapacknvaris:
	cd pfapack/fortran;module purge;module load $(NVARMODULES);cp makefile.nvaris makefile;make clean;make; cp *.mod lib*.a ../../lib/nvaris/; cp *.mod ../../; $(AR) rc ../../lib/nvaris/libMatrix.a *.o; $(AR) rc ../../libMatrix.a *.o

#-----------------------------------------------------------------------------------------------------------------------------
# NVIDIA + gpus on ICARUS:
icarus:   
	make clean
	make pfapackicarus
	module purge;module load $(ICMODULES); make --no-print-directory FC="$(ICFC)" FFLAGS="$(ICFFLAGS)" LIBS="$(ICLIBS)"  mods
	make  libs
	cp libMatrix.a *.mod lib/icarus

pfapackicarus:
	cd pfapack/fortran;module purge;module load $(ICMODULES);cp makefile.icarus makefile;make clean;make; cp *.mod lib*.a ../../lib/icarus/; cp *.mod ../../; $(AR) rc ../../lib/icarus/libMatrix.a *.o; $(AR) rc ../../libMatrix.a *.o

#-----------------------------------------------------------------------------------------------------------------------------
# NVIDIA + gpus on ICARUS:
icarusnogpu:   
	make clean
	make pfapackicarusnogpu
	module purge;module load $(ICPUMODULES); make --no-print-directory FC="$(ICPUFC)" FFLAGS="$(ICPUFFLAGS)" LIBS="$(ICPULIBS)"  mods
	make  libs
	cp libMatrix.a *.mod lib/icarusnogpu

pfapackicarusnogpu:
	cd pfapack/fortran;module purge;module load $(ICPUMODULES);cp makefile.icarusnogpu makefile;make clean;make; cp *.mod lib*.a ../../lib/icarusnogpu/; cp *.mod ../../; $(AR) rc ../../lib/icarusnogpu/libMatrix.a *.o; $(AR) rc ../../libMatrix.a *.o

#-----------------------------------------------------------------------------------------------------------------------------
gnu:    
	make clean
	make pfapack
	make --no-print-directory FC="$(FC)" FFLAGS="$(FFLAGS)" LIBS="$(LIBS)"  mods
	$(AR) rc libMatrix.a $(OBJS)
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
	$(AR) rcv libMatrix.a $(OBJS)

funclist:
	@awk '/:\.:\.:/{print substr($$0,0,138)}' $(SRCS) >  matrix_proclist;a2ps -1 --chars-per-line=138 -r -o - matrix_proclist | ps2pdf - > matrix_proclist.pdf

clean:
	$(RM) test test_benchmark test.dat *.o *.mod core libMatrix.a fort.[0-9][0-9] Matrix.tgz libpfapack.a

pack:
	@./pack
