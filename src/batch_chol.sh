#PBS -N Cholesky
#PBS -l nodes=1:towel:ppn=12
#PBS -o pbs_output.out
#PBS -e pbs_error.out
### send me email when job begins
#PBS -m b
### send me email when job end
#PBS -m e
### send me email when job aborts (with an error)
#PBS -m a

D=100000
NP=12

cd /nfs/user01/jimmie21/github/hfrisk/src
 
LIBDIR=/nfs/user01/jimmie21/github/hfrisk/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${LIBDIR}/gsl-1.15/install/gsl/lib:${LIBDIR}/nlopt-2.3/install/lib:${LIBDIR}/armadillo-3.6.1/install/usr/lib64:${LIBDIR}/fftw-3.3.2/installdir/lib:${LIBDIR}/hdf5-1.8.10/install/lib

OUTFILE=/nfs/user01/jimmie21/github/hfrisk/src/runTimes/chol_d=${D}_np=${NP}.txt

mpirun -np ${NP} ../bin/Cholesky -n ${D} > ${OUTFILE}

