
CG_DATA="/Users/jimmiegoode/Documents/Glimm/github/hfrisk/data/csi/patentData"
CG_LUT="/Users/jimmiegoode/Documents/Glimm/github/hfrisk/data/LUTs"
CG_EXP="${CG_DATA}/exports"

mpirun -np 2 ../bin/HFRisk \
	-df "${CG_DATA}/csi_20030101_20120801_v3" \
	-rf "${CG_EXP}/nill4" \
	-mc 5 \
	-nSims 10000 \
	-v 1 \
	-margOnly 0 \
	-innovType 2 \
	-doCheckEigs 0 \
	-beginDate 20080102 \
	-endDate 20080105 \
	-lut_path "${CG_LUT}/test1" \
	-doLUT 1 \
	-goBig 0
