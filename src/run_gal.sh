
CG_DATA="/nfs/user03/jimmie21/data"
CG_LUT="/nfs/user03/jimmie21/LUTs"
CG_EXP="${CG_DATA}/exports"

mpirun -np 33 ../bin/HFRisk \
	-x OMP_NUM_THREADS \
	-df "${CG_DATA}/csi_20030101_20120801_v3" \
	-rf "${CG_EXP}/nill" \
	-mc -1 \
	-nSims 10000 \
	-v 0 \
	-margOnly 0 \
	-innovType 2 \
	-depStruct 1 \
	-doCheckEigs 0 \
	-beginDate 20080102 \
	-endDate 20080105 \
	-lut_path "${CG_LUT}/test2" \
	-doLUT 1 \
	-margOnly 1 \
	-goBig 0

