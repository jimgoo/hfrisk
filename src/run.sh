

CG_DATA="/Users/jimmiegoode/Documents/Glimm/github/hfrisk/data/csi/patentData"
CG_LUT="/Users/jimmiegoode/Documents/Glimm/github/hfrisk/data/LUTs"

################################
CG_EXP="${CG_DATA}/exports"

mkdir ${CG_EXP}

mpirun -np 2 ../lib/tester \
	-x OMP_NUM_THREADS \
	-df "${CG_DATA}/csi_20030101_20120801_v3" \
	-rf "${CG_EXP}/20130105_ST" \
	-iS 5 \
	-nSims 10000 \
	-v 0 \
	-margOnly 0 \
	-innovType 1 \
	-depStruct 1 \
	-doCheckEigs 0 \
	-beginDate 20080102 \
	-endDate 20080104 \
	-lut_path "${CG_LUT}/test1" \
	-doLUT 1
