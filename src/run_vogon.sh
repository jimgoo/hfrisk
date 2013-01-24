
LIBDIR=/nfs/user01/jimmie21/github/hfrisk/lib

export LD_LIBRARY_PATH=${LIBDIR}/gsl-1.15/install/gsl/lib:${LIBDIR}/nlopt-2.3/install/lib:${LIBDIR}/armadillo-3.6.1/install/usr/lib64:${LIBDIR}/fftw-3.3.2/installdir/l\
ib:${LIBDIR}/hdf5-1.8.10/install/lib

CG_DATA="/nfs/user01/jimmie21/github/hfrisk/data"
CG_LUT="/nfs/user01/jimmie21/LUTs"
CG_EXP="${CG_DATA}/exports"

# number of processors
NP=252

MC=5109
IT=1

mpirun -np "${NP}" ../bin/HFRisk \
	-df "${CG_DATA}/csi_20030101_20120801_v3" \
	-rf "${CG_EXP}/np=${NP}_mc=${MC}_innov=${IT}" \
	-mc -1 \
	-v 0 \
	-margOnly 1 \
	-innovType "${IT}" \
	-beginDate 20080102 \
	-endDate 20080105 \
	-goBig 0

IT=2

mpirun -np "${NP}" ../bin/HFRisk \
	-df "${CG_DATA}/csi_20030101_20120801_v3" \
	-rf "${CG_EXP}/np=${NP}_mc=${MC}_innov=${IT}" \
	-mc -1 \
	-v 0 \
	-margOnly 1 \
	-innovType "${IT}" \
	-beginDate 20080102 \
	-endDate 20080105 \
	-goBig 0

######################

MC=100000
IT=1

mpirun -np "${NP}" ../bin/HFRisk \
	-df "${CG_DATA}/csi_20030101_20120801_v3" \
	-rf "${CG_EXP}/np=${NP}_mc=${MC}_innov=${IT}" \
	-mc -1 \
	-v 0 \
	-margOnly 1 \
	-innovType "${IT}" \
	-beginDate 20080102 \
	-endDate 20080105 \
	-goBig 100000

IT=2

mpirun -np "${NP}" ../bin/HFRisk \
	-df "${CG_DATA}/csi_20030101_20120801_v3" \
	-rf "${CG_EXP}/np=${NP}_mc=${MC}_innov=${IT}" \
	-mc -1 \
	-v 0 \
	-margOnly 1 \
	-innovType "${IT}" \
	-beginDate 20080102 \
	-endDate 20080105 \
	-goBig 100000

