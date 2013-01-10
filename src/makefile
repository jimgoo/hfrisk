#-------------------------------------------------------------------------------
##  Author: Jimmie Goode
##  Created: 2012-09-01
##  Version: 1.0
#
# mpirun -np 2 ./tester -df ../data/csi_20030101_20120801_v3 -mc 5 -maxT 1009 -nSims 10000 -v 0 -rf ../data/exports/test -margOnly 0 -innovType 3 -depStruct 1
#-------------------------------------------------------------------------------

OBJS     = Constants.o gamma.o IO.o Stats.o skeweduvstudenttpdf.o myskewt.o mygarch.o msj.o mystable.o myskewt_lut.o
#OBJS     = Stats.o mystable.o myassg.o #Constants.o gamma.o IO.o Stats.o skeweduvstudenttpdf.o myskewt.o
#OBJS     = Constants.o Stats.o myskewt.o skeweduvstudenttpdf.o IO.o gamma.o myskewt_lut.o
EXE      = tester
DEBUG    =
CFLAGS   = -O2 -w -c $(DEBUG)
LFLAGS   = -Wall $(DEBUG)
OUTPUT   = -o $(EXE)
LNK      = -lgsl -lgslcblas -lnlopt -lfftw3 -lm

# Link Flag: -lhdf5
# HDF5_I = -I/Users/jimmiegoode/Documents/Glimm/lib/hdf5-1.8.9/install/include
# HDF5_L = -L/Users/jimmiegoode/Documents/Glimm/lib/hdf5-1.8.9/install/lib
HDF5_I = 
HDF5_L = 

#-------------------------------------------------------------------------------
# Target Speciffic Parameters
#-------------------------------------------------------------------------------

TARGET = MAC

#Mac speciffic
ifeq ($(TARGET), MAC)
CC    = mpic++ 
BOOST = -I/opt/local/include
GSL   = -I/opt/local/include/gsl
MPI   = -I/opt/local/include/openmpi
NLOPT = -I/usr/local/include
AD_I  = -I/Users/jimmiegoode/Documents/Glimm/lib/armadillo-3.4.2/include
AD_L  =  -framework Accelerate
FFTW_I = -I/Users/jimmiegoode/Documents/Glimm/lib/fftw-3.3.2/installdir/include
FFTW_L = -L/Users/jimmiegoode/Documents/Glimm/lib/fftw-3.3.2/installdir/lib/
endif

# Galaxy speciffic
ifeq ($(TARGET), GALAXY)
CC    = mpic++  -fopenmp
BOOST = -I /nfs/user03/copula/20120323/lib/boost_1_49_0
GSL   = -I /nfs/user03/copula/20120323/lib/gsl-1.15/INSTALL_DIR/gsl/include/gsl
MPI   = -I /usr/local/pkg/openmpi/include
NLOPT = -I/nfs/user03/copula/20120323/lib/nlopt-2.2.4/build/include -L/nfs/user03/copula/20120323/lib/nlopt-2.2.4/build/lib

AD_I  = -I/nfs/user03/copula/20120323/lib/armadillo-3.4.1/include -I/nfs/user03/copula/20120323/lib/OpenBLAS-v0.2.3-0/xianyi-OpenBLAS-48f075c/install/include
AD_L  = -L/nfs/user03/copula/20120323/lib/OpenBLAS-v0.2.3-0/xianyi-OpenBLAS-48f075c/install/lib -lopenblas /usr/lib/liblapack.so -lgfortran

FFTW_I = -I/nfs/user03/copula/20120323/lib/fftw-3.3.2/installdir/include
FFTW_L = -I/nfs/user03/copula/20120323/lib/fftw-3.3.2/installdir/lib

# Elemental
include /nfs/user03/jimmie21/libs/conf/elemvariables

#SVD: SVD.cpp
#    ${CXX} ${ELEM_COMPILE_FLAGS} $< -o $@ ${ELEM_LINK_FLAGS} ${ELEM_LIBS}
# ELEM_I = -I/nfs/user03/jimmie21/libs/include
# ELEM_L = -L/nfs/user03/jimmie21/libs/lib

endif


ifeq ($(TARGET), BLUEGENE)
CC    = mpicxx  -fopenmp
BOOST = -I/gpfs/home3/j/jgoode/lib/boost_1_51_0
GSL   = -I/gpfs/home3/j/jgoode/lib/gsl-1.15/install/gsl/include -L/gpfs/home3/j/jgoode/lib/gsl-1.15/install/gsl/lib
MPI   = -I/bgsys/drivers/ppcfloor/comm/include
NLOPT = -I/gpfs/home3/j/jgoode/lib/nlopt-2.3/install/include -L/gpfs/home3/j/jgoode/lib/nlopt-2.3/install/lib

#AD = -I/gpfs/home3/j/jgoode/lib/armadillo-3.4.2/include -L/gpfs/home3/j/jgoode/examples/lapack/lapack-3.3.0

AD_I = -I/gpfs/home3/j/jgoode/lib/armadillo-3.4.2/include
AD_L = -L/gpfs/home3/j/jgoode/examples/lapack/lapack-3.3.0
#AD_FLAG  = -framework Accelerate #-larmadillo  #-llapack -lblas -lgfortran
endif


#export OMP_NUM_THREADS=2
#echo "OMP_NUM_THREADS = $OMP_NUM_THREADS"


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

all: main

main: $(OBJS)
	$(CC) $(LFLAGS) $(OUTPUT) $(OBJS) $(GSL) $(BOOST) $(NLOPT) $(AD_L) $(MPI) $(FFTW_L) $(HDF5_L) $(LNK)

Constants.o: Constants.h Constants.cpp
	$(CC) $(CFLAGS) Constants.cpp

myskewt.o: myskewt.hpp myskewt.cpp Constants.h gamma.h IO.h Stats.hpp
	$(CC) $(CFLAGS) $(BOOST) $(GSL) $(NLOPT) $(AD_I) myskewt.cpp

msj.o: msj.h msj.cpp IO.h Stats.hpp mygarch.hpp myskewt.hpp mystable.hpp myskewt_lut.hpp
	$(CC) $(CFLAGS) $(GSL) $(BOOST) $(MPI) $(NLOPT) $(AD_I) msj.cpp

IO.o: IO.h IO.cpp
	$(CC) $(CFLAGS) $(GSL) $(BOOST) IO.cpp

Stats.o: Stats.hpp Stats.cpp
	$(CC) $(CFLAGS) $(GSL) $(BOOST) $(AD_I) $(MPI) $(HDF5_I) Stats.cpp

gamma.o: gamma.h gamma.cpp
	$(CC) $(CFLAGS) $(GSL) $(BOOST) gamma.cpp

skeweduvstudenttpdf.o: skeweduvstudenttpdf.h skeweduvstudenttpdf.cpp gamma.h
	$(CC) $(CFLAGS) $(GSL) $(BOOST) skeweduvstudenttpdf.cpp

mygarch.o: mygarch.hpp mygarch.cpp Stats.hpp IO.h 
	$(CC) $(CFLAGS) $(GSL) $(BOOST) $(AD_I) $(NLOPT) mygarch.cpp

mystable.o: mystable.hpp mystable.cpp Stats.hpp
	$(CC) $(CFLAGS) $(GSL) $(BOOST) $(AD_I) $(NLOPT) $(FFTW_I) mystable.cpp

myassg.o: myassg.cpp Stats.hpp mystable.hpp 
	$(CC) $(CFLAGS) $(GSL) $(BOOST) $(AD_I) myassg.cpp

myskewt_lut.o: myskewt_lut.hpp myskewt_lut.cpp myskewt.hpp
	$(CC) $(CFLAGS) $(GSL) $(BOOST) $(AD_I) $(NLOPT) myskewt_lut.cpp

clean:
	rm -rf *.o *.out $(EXE)

cleantmp:
	rm -f *~

pack:
	rm -f src.zip
	zip -rq src.zip *

send:
	curl -T src.zip ftp://hjkm.net/private/ --user jimmie1921:Geaux235711!

#scp src.zip jgoode@ssh.ams.sunysb.edu:~/dropbox/
#ssh jgoode@ssh.ams.sunysb.edu "sh pushbox"
#scp -i /home/fermat/grads/jgoode/.ssh/galaxy2 -r dropbox/* jimmie21@galaxy.ams.sunysb.edu:~/dropbox/

sltest:
	$(CC) -O2 $(LFLAGS) $(BLAS) $(LAPACK) $(SCALAPACK) $(MPI) sltest.cpp -o sltest 
