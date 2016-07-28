#!/usr/bin/env python

import sys, os

from seqcluster.libs.fastq import collapse, splitext_plus, write_output
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.utils import (file_exists, append_stem, replace_directory, symlink_plus, local_path_export)

from collections import Counter


if __name__ == "__main__":
	in_file = sys.argv[1]

	"""
	Collpase reads into unique sequences with seqcluster
	"""
	out_file = append_stem(in_file, ".trimming").replace(".gz", "")
	#out_file = splitext_plus(os.path.basename(fastq))[0] + ".fq"
	seqs = collapse(in_file)
	write_output(out_file, seqs, 1)

	"""
	Calculate size distribution after adapter removal
	"""
	data = Counter()
	out_stat_file = out_file + "_size_stats"

	with open(out_file) as in_handle:
		for line in in_handle:
			counts = int(line.strip().split("_x")[1])
			line = in_handle.next()
			l = len(line.strip())
			in_handle.next()
			in_handle.next()
			data[l] += counts
	with file_transaction(out_stat_file) as tx_out_stat_file:
		with open(tx_out_stat_file, 'w') as out_handle:
			for l, c in data.iteritems():
				out_handle.write("%s %s\n" % (l, c))	
