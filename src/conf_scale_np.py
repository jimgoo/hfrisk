
import os, codecs
#import numpy as np

#===================================================================


sToday    = "20121018"
innovTypes = [1,2,3,4] #"3"

#nps = [2, 4, 8, 16, 32, 64, 128]
nps = [128]
print "Number of processes to use:\n"
print(nps)


sFile = "conf_run.sh"
fout = codecs.open(sFile, encoding='utf-8', mode='w')
fout.write("# This script was created by 'conf_runner.py'\n\n");


cmd = "export OMP_NUM_THREADS=2\n\n"
fout.write(cmd)
print(cmd)


for innovType in innovTypes:

    reportDir = "../data/exports/" + sToday + "/innovType=" + str(innovType)

    if os.path.exists(reportDir) == False:
        os.makedirs(reportDir)

    for np in nps:
        
        rf = reportDir + "/np=" + str(np)

        cmd = "mpirun -np " + str(np) + " -x OMP_NUM_THREADS ./tester" \
            + " -df ../data/csi_20030101_20120801_v3" \
            + " -mc -1" \
            + " -maxT 1009" \
            + " -nSims 10000" \
            + " -v 0" \
            + " -rf " + rf \
            + " -margOnly 1" \
            + " -innovType " + str(innovType) \
            + "\n"

        fout.write(cmd)
        print(cmd)
	
fout.close()
print "Script saved as " + sFile + "\n"
