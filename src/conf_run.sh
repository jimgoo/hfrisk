# This script was created by 'conf_runner.py'

export OMP_NUM_THREADS=2

mpirun -np 128 -x OMP_NUM_THREADS ./tester -df ../data/csi_20030101_20120801_v3 -mc 1000 -maxT 1009 -nSims 10000 -v 0 -rf ../data/exports/20121020/innovType=1/n=1000 -margOnly 0 -innovType 1
mpirun -np 128 -x OMP_NUM_THREADS ./tester -df ../data/csi_20030101_20120801_v3 -mc 2000 -maxT 1009 -nSims 10000 -v 0 -rf ../data/exports/20121020/innovType=1/n=2000 -margOnly 0 -innovType 1
mpirun -np 128 -x OMP_NUM_THREADS ./tester -df ../data/csi_20030101_20120801_v3 -mc 3000 -maxT 1009 -nSims 10000 -v 0 -rf ../data/exports/20121020/innovType=1/n=3000 -margOnly 0 -innovType 1
mpirun -np 128 -x OMP_NUM_THREADS ./tester -df ../data/csi_20030101_20120801_v3 -mc 4000 -maxT 1009 -nSims 10000 -v 0 -rf ../data/exports/20121020/innovType=1/n=4000 -margOnly 0 -innovType 1
mpirun -np 128 -x OMP_NUM_THREADS ./tester -df ../data/csi_20030101_20120801_v3 -mc 5000 -maxT 1009 -nSims 10000 -v 0 -rf ../data/exports/20121020/innovType=1/n=5000 -margOnly 0 -innovType 1
