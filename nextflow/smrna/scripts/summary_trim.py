import os, sys

from collections import Counter

from bcbio.utils import (file_exists, append_stem, replace_directory, symlink_plus, local_path_export)

if __name__ == "__main__":
	"""
	Calculate size distribution after adapter removal
	"""
	in_file = '/BiO/BioPeople/brandon/test_nextflow/smrna/outdir/trimmed/SRR950876.fastq/SRR950876_trimmed.fq'

	data = Counter()
	out_file = in_file + "_size_stats"
	if file_exists(out_file):
		sys.exit()

	with open(in_file) as in_handle:
		for line in in_handle:
			counts = int(line.strip().split("_x")[1])
			line = in_handle.next()
			l = len(line.strip())
			in_handle.next()
			in_handle.next()
			data[l] += counts
	with file_transaction(out_file) as tx_out_file:
		with open(tx_out_file, 'w') as out_handle:
			for l, c in data.iteritems():
				out_handle.write("%s %s\n" % (l, c))
