#-------------------------------------------------------------------------------
##  Author: Jimmie Goode
##  Created: 2012-09-01
#-------------------------------------------------------------------------------

OBJS = Constants.o gamma.o IO.o Stats.o skeweduvstudenttpdf.o myskewt.o mygarch.o msj.o mystable.o myskewt_lut.o
#OBJS = Stats.o mystable.o myassg.o #Constants.o gamma.o IO.o Stats.o skeweduvstudenttpdf.o myskewt.o
#OBJS = Constants.o Stats.o myskewt.o skeweduvstudenttpdf.o IO.o gamma.o myskewt_lut.o
#OBJS     = mygarch_par.o #MasterSlave.o
#OBJS     = Backtest.o MasterSlave.o mygarch_par.o IO.o
#OBJS     = Copula.o Stats.o

#OBJS = Stats.o myskewt.o IO.o skeweduvstudenttpdf.o gamma.o mygarch.o

BINDIR   = ../bin
EXE      = $(BINDIR)/HFRisk
DEBUG    =
CFLAGS   = -O2 -w -c $(DEBUG)
LFLAGS   = -Wall $(DEBUG)
LNK      = -lgsl -lgslcblas -lnlopt -lfftw3 -lm

# extra linker flags for MAC mpic++ only.
MPICCF = 

#-------------------------------------------------------------------------------
# Target Speciffic Parameters
#-------------------------------------------------------------------------------

#TARGET = MAC
# TARGET = {MAC, GALAXY, VOGON}
#
# Pass this with
#    $ make TARGET=MAC

#Mac speciffic
ifeq ($(TARGET), MAC)
#CCMPI = mpic++ # uses clang++ which is broke
CCMPI  = g++
MPICCF = -L/opt/local/lib -lmpi_cxx -lmpi -lm
BOOST  = -I/opt/local/include
GSL    = -I/opt/local/include/gsl
MPI    = -I/opt/local/include/openmpi
NLOPT  = -I/usr/local/include
AD_I   = -I/Users/jimmiegoode/Documents/Glimm/lib/armadillo-3.6.1/include
AD_L   = -framework Accelerate
FFTW_I = -I/Users/jimmiegoode/Documents/Glimm/lib/fftw-3.3.2/installdir/include
FFTW_L = -L/Users/jimmiegoode/Documents/Glimm/lib/fftw-3.3.2/installdir/lib/
HDF_I  = -I/Users/jimmiegoode/Documents/Glimm/lib/hdf5-1.8.9/install/include
HDF_L  = -L/Users/jimmiegoode/Documents/Glimm/lib/hdf5-1.8.9/install/lib -lhdf5

include /Users/jimmiegoode/Documents/Glimm/lib/elemental/elemental-0.77/build/conf/elemvariables

endif

# Galaxy speciffic
ifeq ($(TARGET), GALAXY)
CCMPI    = mpic++  -fopenmp
BOOST = -I /nfs/user03/copula/20120323/lib/boost_1_49_0
GSL   = -I /nfs/user03/copula/20120323/lib/gsl-1.15/INSTALL_DIR/gsl/include/gsl
MPI   = -I /usr/local/pkg/openmpi/include
NLOPT = -I/nfs/user03/copula/20120323/lib/nlopt-2.2.4/build/include -L/nfs/user03/copula/20120323/lib/nlopt-2.2.4/build/lib

# Use OpenBLAS for threading and regular built in LAPACK (/usr/lib/liblapack.so)
AD_I  = -I/nfs/user03/copula/20120323/lib/armadillo-3.6.1/include -I/nfs/user03/copula/20120323/lib/OpenBLAS-v0.2.3-0/xianyi-OpenBLAS-48f075c/install/include
AD_L  = -L/nfs/user03/copula/20120323/lib/OpenBLAS-v0.2.3-0/xianyi-OpenBLAS-48f075c/install/lib -lopenblas /usr/lib/liblapack.so -lgfortran
FFTW_I = -I/nfs/user03/copula/20120323/lib/fftw-3.3.2/installdir/include
FFTW_L = -I/nfs/user03/copula/20120323/lib/fftw-3.3.2/installdir/lib
HDF_I = -I/nfs/user03/copula/20120323/lib/hdf5-1.8.9/install/include
HDF_L = -L/nfs/user03/copula/20120323/lib/hdf5-1.8.9/install/lib -lhdf5

include /nfs/user03/copula/20120323/lib/elemental-0.77/build/conf/elemvariables

endif


ifeq ($(TARGET), BLUEGENE)
CCMPI    = mpicxx  -fopenmp
BOOST = -I/gpfs/home3/j/jgoode/lib/boost_1_51_0
GSL   = -I/gpfs/home3/j/jgoode/lib/gsl-1.15/install/gsl/include -L/gpfs/home3/j/jgoode/lib/gsl-1.15/install/gsl/lib
MPI   = -I/bgsys/drivers/ppcfloor/comm/include
NLOPT = -I/gpfs/home3/j/jgoode/lib/nlopt-2.3/install/include -L/gpfs/home3/j/jgoode/lib/nlopt-2.3/install/lib

#AD = -I/gpfs/home3/j/jgoode/lib/armadillo-3.4.2/include -L/gpfs/home3/j/jgoode/examples/lapack/lapack-3.3.0

AD_I = -I/gpfs/home3/j/jgoode/lib/armadillo-3.4.2/include
AD_L = -L/gpfs/home3/j/jgoode/examples/lapack/lapack-3.3.0
#AD_FLAG  = -framework Accelerate #-larmadillo  #-llapack -lblas -lgfortran
endif


ifeq ($(TARGET), VOGON)

LIBDIR=/nfs/user01/jimmie21/github/hfrisk/lib

# ALL libs must be added to the LD_LIBRARY_PATH for Vogon

#export LD_LIBRARY_PATH=${LIBDIR}/gsl-1.15/install/gsl/lib:${LIBDIR}/nlopt-2.3/install/lib:${LIBDIR}/armadillo-3.6.1/install/usr/lib64:${LIBDIR}/fftw-3.3.2/installdir/lib:${LIBDIR}/hdf5-1.8.10/install/lib

CCMPI  = mpic++
BOOST  = -I${LIBDIR}/boost_1_52_0
GSL    = -I${LIBDIR}/gsl-1.15/install/gsl/include -L${LIBDIR}/gsl-1.15/install/gsl/lib
MPI    = -I/usr/local/pkg/openmpi/include
NLOPT  = -I${LIBDIR}/nlopt-2.3/install/include -L${LIBDIR}/nlopt-2.3/install/lib

AD_I   = -I${LIBDIR}/armadillo-3.6.1/install/usr/include
AD_L   = -L${LIBDIR}/armadillo-3.6.1/install/usr/lib64 -larmadillo

FFTW_I = -I${LIBDIR}/fftw-3.3.2/installdir/include
FFTW_L = -L${LIBDIR}/fftw-3.3.2/installdir/lib

HDF_I  = -I${LIBDIR}/hdf5-1.8.10/install/include
HDF_L  = -L${LIBDIR}/hdf5-1.8.10/install/lib -lhdf5

include ${LIBDIR}/elemental-0.77/install/conf/elemvariables


endif

#export OMP_NUM_THREADS=2
#echo "OMP_NUM_THREADS = $OMP_NUM_THREADS"


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

all: main

main: $(OBJS)
	$(CCMPI) $(LFLAGS) -o $(EXE) $(OBJS) $(GSL) $(BOOST) $(NLOPT) $(HDF_L) $(AD_L) $(MPI) $(FFTW_L) $(LNK) ${ELEM_LINK_FLAGS} ${ELEM_LIBS} $(MPICCF) ${ELEM_LINK_FLAGS} ${ELEM_LIBS}

Constants.o: Constants.h Constants.cpp
	$(CCMPI) $(CFLAGS) Constants.cpp

myskewt.o: myskewt.hpp myskewt.cpp Constants.h gamma.h IO.hpp Stats.hpp
	$(CCMPI) $(CFLAGS) $(GSL) $(BOOST) $(HDF_I) $(MPI) $(NLOPT) $(AD_I) myskewt.cpp

msj.o: msj.h msj.cpp IO.hpp Stats.hpp mygarch.hpp myskewt.hpp mystable.hpp myskewt_lut.hpp
	$(CCMPI) $(CFLAGS) $(GSL) $(BOOST) $(HDF_I) $(MPI) $(NLOPT) $(AD_I) $(FFTW_I) msj.cpp

IO.o: IO.hpp IO.cpp
	$(CCMPI) $(CFLAGS) $(GSL) $(BOOST) $(HDF_I) $(AD_I) IO.cpp

Stats.o: Stats.hpp Stats.cpp
	$(CCMPI) $(CFLAGS) $(GSL) $(BOOST) $(HDF_I) $(AD_I) $(MPI) $(HDF5_I) Stats.cpp

gamma.o: gamma.h gamma.cpp
	$(CCMPI) $(CFLAGS) $(GSL) $(BOOST) gamma.cpp

skeweduvstudenttpdf.o: skeweduvstudenttpdf.h skeweduvstudenttpdf.cpp gamma.h
	$(CCMPI) $(CFLAGS) $(GSL) $(BOOST) skeweduvstudenttpdf.cpp

mygarch.o: mygarch.hpp mygarch.cpp Stats.hpp IO.hpp 
	$(CCMPI) $(CFLAGS) $(GSL) $(BOOST) $(HDF_I) $(MPI) $(AD_I) $(NLOPT) mygarch.cpp

mystable.o: mystable.hpp mystable.cpp Stats.hpp
	$(CCMPI) $(CFLAGS) $(GSL) $(BOOST) $(HDF_I) $(MPI) $(AD_I) $(NLOPT) $(FFTW_I) mystable.cpp

myassg.o: myassg.cpp Stats.hpp mystable.hpp 
	$(CCMPI) $(CFLAGS) $(GSL) $(BOOST) $(HDF_I) $(MPI) $(AD_I) myassg.cpp

myskewt_lut.o: myskewt_lut.hpp myskewt_lut.cpp myskewt.hpp
	$(CCMPI) $(CFLAGS) $(GSL) $(BOOST) $(HDF_I) $(MPI) $(AD_I) $(NLOPT) myskewt_lut.cpp

ShardBinaries: ShardBinaries.cpp
	g++ -o $(BINDIR)/ShardBinaries ShardBinaries.cpp $(BOOST) $(HDF_I) $(HDF_L) $(AD_I) $(AD_L)

Csv2Binary: Csv2Binary.cpp
	g++ -o $(BINDIR)/Csv2Binary Csv2Binary.cpp $(BOOST) $(HDF_I) $(HDF_L) $(AD_I) $(AD_L)

ArmaTester: ArmaTester.cpp
	g++ -o $(BINDIR)/ArmaTester ArmaTester.cpp $(AD_I) $(AD_L)

mygarch_par.o: mygarch_par.cpp mygarch.hpp MasterSlave.hpp
	$(CCMPI) $(CFLAGS) ${ELEM_COMPILE_FLAGS} $(MPI) $(AD_I) mygarch_par.cpp

MasterSlave.o: MasterSlave.hpp MasterSlave.cpp
	$(CCMPI) $(CFLAGS) $(MPI) MasterSlave.cpp

SymmetricEig: SymmetricEig.cpp
	${CXX} ${ELEM_COMPILE_FLAGS} $< -o $@ ${ELEM_LINK_FLAGS} ${ELEM_LIBS}

Backtest.o: Backtest.hpp Backtest.cpp
	$(CCMPI) $(CFLAGS) $(MPI) $(AD_I) ${ELEM_COMPILE_FLAGS} Backtest.cpp

Copula.o: Copula.cpp Stats.o
	$(CCMPI) $(CFLAGS) $(GSL) $(BOOST) $(HDF_I) $(MPI) $(NLOPT) $(AD_I) ${ELEM_COMPILE_FLAGS} Copula.cpp


# Dependency Estimation
dep: Copula.o Stats.o
	$(CCMPI) $(LFLAGS) -o ../bin/DepStruct Copula.o Stats.o $(GSL) $(BOOST) $(NLOPT) $(HDF_L) $(AD_L) $(MPI) $(FFTW_L) $(LNK) ${ELEM_LINK_FLAGS} ${ELEM_LIBS} $(MPICCF) ${ELEM_LINK_FLAGS} ${ELEM_LIBS}


# Cholesky Factorization
Cholesky: Cholesky.cpp
	${CXX} ${ELEM_COMPILE_FLAGS} $(HDF_I) $(HDF_L) $(AD_I) $(AD_L) $< -o ../bin/$@ ${ELEM_LINK_FLAGS} ${ELEM_LIBS}


# Cholesky Factorization - Elemental
CholeskyElem: CholeskyElem.cpp
	${CXX} ${ELEM_COMPILE_FLAGS} $< -o ../bin/$@ ${ELEM_LINK_FLAGS} ${ELEM_LIBS}


# Monte Carlo VaR
MC: MonteCarlo.cpp
	${CXX} ${ELEM_COMPILE_FLAGS} $(AD_I) $(AD_L) $< -o ../bin/$@ ${ELEM_LINK_FLAGS} ${ELEM_LIBS} $(BOOST)


clean:
	rm -rf *.o *.out $(EXE)

cleantmp:
	rm -f *~

pack:
	rm -f src.zip
	zip -rq src.zip *

send:
	curl -T src.zip ftp://hjkm.net/private/ --user jimmie1921

sltest:
	$(CCMPI) -O2 $(LFLAGS) $(BLAS) $(LAPACK) $(SCALAPACK) $(MPI) sltest.cpp -o sltest 

test:
	echo $(TARGET)
