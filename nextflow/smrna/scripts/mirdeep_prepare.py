import os
import sys
import os.path as op

import pysam

sys.path.append('/BiO/BioTools/bcbio/data/anaconda/lib/python2.7/site-packages/bcbio')
from bcbio.log import logger
from bcbio.utils import file_exists, safe_makedir, chdir, get_perl_exports
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils


def _prepare_inputs(ma_fn, bam_file, out_dir):
    """
    Convert to fastq with counts
    """
    fixed_fa = os.path.join(out_dir, "file_reads.fa")
    count_name =dict()
    with file_transaction(fixed_fa) as out_tx:
        with open(out_tx, 'w') as out_handle:
            with open(ma_fn) as in_handle:
                h = in_handle.next()
                for line in in_handle:
                    cols = line.split("\t")
                    name_with_counts = "%s_x%s" % (cols[0], sum(map(int, cols[2:])))
                    count_name[cols[0]] = name_with_counts
                    print >>out_handle, ">%s\n%s" % (name_with_counts, cols[1])
    fixed_bam = os.path.join(out_dir, "align.bam")
    bam_handle = pysam.AlignmentFile(bam_file, "rb")
    with pysam.AlignmentFile(fixed_bam, "wb", template=bam_handle) as out_handle:
        for read in bam_handle.fetch():
            read.query_name = count_name[read.query_name]
            out_handle.write(read)

    return fixed_fa, fixed_bam


if __name__ == "__main__":
	bam_file = sys.argv[1] # seqs.bam
	collapsed = sys.argv[2] # seqs.ma

	out_dir = "./"

	collapsed, bam_file = _prepare_inputs(collapsed, bam_file, out_dir)


