OUTFILE=/nfs/user01/jimmie21/github/hfrisk/src/RT_DEP
mkdir ${OUTFILE}

LIBDIR=/nfs/user01/jimmie21/github/hfrisk/lib

export LD_LIBRARY_PATH=${LIBDIR}/gsl-1.15/install/gsl/lib:${LIBDIR}/nlopt-2.3/install/lib:${LIBDIR}/armadillo-3.6.1/install/usr/lib64:${LIBDIR}/fftw-3.3.2/installdir/l\
ib:${LIBDIR}/hdf5-1.8.10/install/lib

mpirun -np 204 ../bin/Cholesky -n 100000
