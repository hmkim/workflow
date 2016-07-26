#!/usr/bin/env python
#https://github.com/chapmanb/bcbio-nextgen/blob/4c2480c32bc6833b0ade41aa518b9b650196c3bc/bcbio/srna/group.py


import os
import os.path as op
import sys
import shutil
from collections import namedtuple

try:
	from seqcluster import prepare_data as prepare
	from seqcluster import templates as template_seqcluster
	from seqcluster.seqbuster import _create_counts, _read_miraligner, _tab_output
except ImportError:
	pass

from bcbio.utils import file_exists, safe_makedir, move_safe, append_stem
from bcbio.distributed.transaction import file_transaction

ids = sys.argv[1]
files = sys.argv[2]


out_dir = "./"

fn = []

args = namedtuple('args', 'debug print_debug minc minl maxl out')
args = args(False, False, 2, 17, 40, out_dir)

ma_out = op.join(out_dir, "seqs.ma")
seq_out = op.join(out_dir, "seqs.fastq")
min_shared = max(int(len(fn) / 10.0), 1)


arr_ids  = ids.split(',')
arr_files = files.split(',')

for id_, file_ in zip( arr_ids, arr_files ):
	fn.append("%s\t%s"%(file_,id_))


if not file_exists(ma_out):
	seq_l, sample_l = prepare._read_fastq_files(fn, args)
	with file_transaction(ma_out) as ma_tx:
		with open(ma_tx, 'w') as ma_handle:
			with open(seq_out, 'w') as seq_handle:
				prepare._create_matrix_uniq_seq(sample_l, seq_l, ma_handle, seq_handle, min_shared)


