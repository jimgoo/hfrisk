
CG_DATA="/nfs/user03/jimmie21/data"
CG_LUT="/nfs/user03/jimmie21/LUTs"

CG_EXP="${CG_DATA}/exports"

mkdir ${CG_EXP}

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
	-goBig 80000

# mpirun -np 4 ../lib/tester \
# 	-x OMP_NUM_THREADS \
# 	-df /nfs/user03/jimmie21/data/csi_20030101_20120801_v3 \
# 	-rf /nfs/user03/jimmie21/data/exports/20130105_ST \
# 	-mc 5 \
# 	-nSims 10000 \
# 	-v 1 \
# 	-margOnly 0 \
# 	-innovType 1 \
# 	-depStruct 1 \
# 	-doCheckEigs 0 \
# 	-beginDate 20080102 \
# 	-endDate 20080104 \
# 	-lut_path /nfs/user03/jimmie21/LUTs/test2 \
# 	-doLUT 1

# mpirun -np 64 ../lib/tester -x OMP_NUM_THREADS -df /nfs/user03/jimmie21/data/csi_20030101_20120801_v3 -rf /nfs/user03/jimmie21/data/exports/20130105_ST -nSims 10000 -v 1 -margOnly 0 -innovType 1 -depStruct 1 -doCheckEigs 0 -beginDate 20080102 -endDate 20080104 -lut_path /nfs/user03/jimmie21/LUTs/t2 -doLUT 1 -iS 5
