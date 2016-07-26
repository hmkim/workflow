"3""Functions to access json with classes"""
from collections import defaultdict
import pandas as pd
import sys
import os.path as op
import os # for os.path
from bcbio.pipeline import datadict as dd


from pprint import pprint # for debug

try:
    from seqcluster import prepare_data as prepare
    from seqcluster import templates as template_seqcluster
    from seqcluster.seqbuster import _create_counts, _read_miraligner, _tab_output
except ImportError:
    pass


from bcbio.pipeline import datadict as dd
from bcbio.utils import file_exists, safe_makedir, move_safe, append_stem

def _make_isomir_counts(data, srna_type="seqbuster", out_dir=None, stem=""):
	"""
	Parse miraligner files to create count matrix.
	"""
	work_dir = dd.get_work_dir(data[0][0])
	if not out_dir:
		out_dir = op.join(work_dir, "mirbase")
	out_novel_isomir = append_stem(op.join(out_dir, "counts.tsv"), stem)
	out_novel_mirna = append_stem(op.join(out_dir, "counts_mirna.tsv"), stem)
	#print "Create %s count data at %s." % (srna_type, out_dir)
	if file_exists(out_novel_mirna):
		return [out_novel_mirna, out_novel_isomir]
	out_dts = []
	for sample in data:
		if sample[0].get(srna_type):
			miraligner_fn = sample[0].get(srna_type)
			print miraligner_fn
			reads = _read_miraligner(miraligner_fn)
#	for sample, miraligner_fn in zip( arr_ids, arr_files ):
			if reads:
				out_file, dt, dt_pre = _tab_output(reads, miraligner_fn + ".back", sample[1])
				out_dts.append(dt)
			else:
				aa =  "WARNING::%s has NOT miRNA annotated.  Check if fasta files is small or species value."%(sample[1])
	if out_dts:
		out_files = _create_counts(out_dts, out_dir)
		out_files = [move_safe(out_files[0], out_novel_isomir), move_safe(out_files[1], out_novel_mirna)]
		return out_files

if __name__ == "__main__":
	ids = sys.argv[1]

	arr_ids  = ids.split(',')

	data = []
	for sample in arr_ids:
		dic1 = { 'seqbuster' : "%s.mirna"%sample, 'seqbuster_novel' :"%s_novel.mirna"%sample}
		#data.append([os.path.dirname(os.path.abspath(sample)),sample])
		data.append([dic1,sample,os.path.dirname(os.path.abspath(sample))])

	work_dir = ".";
	out_mirna = _make_isomir_counts(data, out_dir=op.join(work_dir, "mirbase"))
	out_novel = _make_isomir_counts(data, "seqbuster_novel", op.join(work_dir, "mirdeep2"), "_novel")


	
