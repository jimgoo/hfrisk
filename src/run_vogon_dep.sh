
LIBDIR=/nfs/user01/jimmie21/github/hfrisk/lib

export LD_LIBRARY_PATH=${LIBDIR}/gsl-1.15/install/gsl/lib:${LIBDIR}/nlopt-2.3/install/lib:${LIBDIR}/armadillo-3.6.1/install/usr/lib64:${LIBDIR}/fftw-3.3.2/installdir/l\
ib:${LIBDIR}/hdf5-1.8.10/install/lib

## 5,109 stocks
#mpirun -np 2 ../bin/DepStruct -garchFile /nfs/user01/jimmie21/github/hfrisk/data/exports/ST_5109/mnResults/20080102.h5

## 100,000 stocks
## 4 nodes X 12 proc = 48 cores
mpirun -np 48 ../bin/DepStruct -garchFile /nfs/user01/jimmie21/github/hfrisk/data/exports/ST_100000/mnResults/20080102.h5 