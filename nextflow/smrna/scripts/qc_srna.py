#https://github.com/chapmanb/bcbio-nextgen/blob/c522bfeb327907994cd5aba137211f1e94ecd462/bcbio/qc/srna.py
import pandas as pd
import os

from bcbio import utils
from bcbio.distributed.transaction import file_transaction, tx_tmpdir


def _get_stats_from_miraligner(fn, out_file, name):
	df = pd.read_csv(fn, sep="\t", dtype={"mism": "string",
		"add": "string",
		"t5": "string",
		"t3": "string"},
		na_values=["."])
	dfmirs = df[['mir', 'freq']].groupby(['mir']).count()
	df5 = df.loc[df.t5 != "0", ['mir', 't5']].groupby(['mir']).count()
	df3 = df.loc[df.t3 != "0", ['mir', 't3']].groupby(['mir']).count()
	dfadd = df.loc[df["add"] != "0", ['mir', 'add']].groupby(['mir']).count()
	dfmut = df.loc[df.mism != "0", ['mir', 'mism']].groupby(['mir']).count()
	with file_transaction(out_file) as tx_out:
		with open(tx_out, "w") as out_handle:
			print >>out_handle, "# stats {name}".format(**locals())
			print >>out_handle, ("mirs\t{mirs}\nisomirs\t{isomirs}").format(
			mirs=len(dfmirs.index), isomirs=len(df.index))
			print >>out_handle, ("mirs_mutations\t{muts}\nmirs_additions\t{add}").format(
			muts=len(dfmut.index), add=len(dfadd.index))
			print >>out_handle, ("mirs_5-trimming\t{t5}\nmirs_3-trimming\t{t3}").format(
			t5=len(df5.index), t3=len(df3.index))
			print >>out_handle, ("iso_mutations\t{muts}\niso_additions\t{add}").format(
			muts=sum(dfmut.mism), add=sum(dfadd["add"]))
			print >>out_handle, ("iso_5-trimming\t{t5}\niso_3-trimming\t{t3}").format(
			t5=sum(df5.t5), t3=sum(df3.t3))

if __name__ == "__main__":
	mirbase_fn = sys.argv[1]
	mirdeep_fn = sys.argv[2]
	#mirbase_fn = '/BiO/BioPeople/brandon/test_nextflow/smrna/work/f3/e4a20e3f6a7d5dcef30af43f336cc5/SRR950879.fastq.mirna'
	#mirdeep_fn = '/BiO/BioPeople/brandon/test_nextflow/smrna/work/4e/a2c1d7cb865b25caf18c23cf85c3cc/SRR950880.fastq_novel.mirna'
	if mirbase_fn:
		out_file = "sample_bcbio_mirbase.txt"
		_get_stats_from_miraligner(mirbase_fn, out_file, "seqbuster")

	if mirdeep_fn:
		out_file = "sample_bcbio_mirdeep2.txt" 
		_get_stats_from_miraligner(mirdeep_fn, out_file, "mirdeep2")


