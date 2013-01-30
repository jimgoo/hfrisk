#PBS -N batch_chol3
#PBS -l nodes=4:towel:ppn=12
#PBS -o /dev/null
#PBS -e /dev/null

D=100000

cd /nfs/user01/jimmie21/github/hfrisk/src

LIBDIR=/nfs/user01/jimmie21/github/hfrisk/lib

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${LIBDIR}/gsl-1.15/install/gsl/lib:${LIBDIR}/nlopt-2.3/install/lib:${LIBDIR}/armadillo-3.6.1/install/usr/lib64:${LIBDIR}/fftw-3.3.2/installdir/lib:${LIBDIR}/hdf5-1.8.10/install/lib

OUTDIR=/nfs/user01/jimmie21/github/hfrisk/src/runTimes2

NP=48
OUTFILE="${OUTDIR}/chol_d=${D}_np=${NP}.txt"
mpirun -np ${NP} ../bin/Cholesky -n ${D} > ${OUTFILE}

