#!/usr/bin/env python

import sys, os
import shutil
from collections import namedtuple

from seqcluster.prepare_data import _read_fastq_files, _create_matrix_uniq_seq

filelist = sys.argv[1]

arg = namedtuple('args', 'minl maxl minc out')
args = arg(15, 40, 1, './')
seq_l, list_s = _read_fastq_files(open(filelist), args)

ma_out = open("seqs.ma", 'w')
seq_out = open("seqs.fa", 'w')
_create_matrix_uniq_seq(list_s, seq_l, ma_out, seq_out, 1)
