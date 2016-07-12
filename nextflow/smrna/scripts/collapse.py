#!/usr/bin/env python

import sys, os

from seqcluster.libs.fastq import collapse, splitext_plus, write_output

fastq = sys.argv[1]

seqs = collapse(fastq)
out_file = splitext_plus(os.path.basename(fastq))[0] + "_trimmed.fastq"
write_output(out_file, seqs, 1)
