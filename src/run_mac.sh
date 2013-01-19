
CG_DATA="/Users/jimmiegoode/Documents/Glimm/github/hfrisk/data/csi/patentData"
CG_LUT="/Users/jimmiegoode/Documents/Glimm/github/hfrisk/data/LUTs"
CG_EXP="${CG_DATA}/exports"

mpirun -np 2 ../bin/HFRisk \
	-x OMP_NUM_THREADS \
	-df "${CG_DATA}/csi_20030101_20120801_v3" \
	-rf "${CG_EXP}/nill" \
	-mc -1 \
	-nSims 10000 \
	-v 1 \
	-margOnly 0 \
	-innovType 1 \
	-depStruct 1 \
	-doCheckEigs 0 \
	-beginDate 20080102 \
	-endDate 20080105 \
	-lut_path "${CG_LUT}/test2" \
	-doLUT 0 \
	-margOnly 1 \
	-goBig 0
