import glob
import sys
import os

in_dir = sys.argv[1]
samplefiles = glob.glob(os.path.join(in_dir, "*.fastq"))
R1 = [x for x in samplefiles if "_R1_" in x]
samplenames = [x.split("_R1_")[0] for x in R1]
for x in samplenames:
    print x
