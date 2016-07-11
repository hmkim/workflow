#!/usr/bin/env python

import sys, os

sys.path.append('/BiO/BioTools/bcbio/data/anaconda/lib/python2.7/site-packages/seqcluster')
from libs.fastq import collapse, splitext_plus, write_output

fastq = sys.argv[1]

seqs = collapse(fastq)
out_file = splitext_plus(os.path.basename(fastq))[0] + "_trimmed.fastq"
write_output(out_file, seqs, 1)
