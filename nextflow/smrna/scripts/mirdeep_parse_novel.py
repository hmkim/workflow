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

def _parse_novel(csv_file, sps="new"):
    """Create input of novel miRNAs from miRDeep2"""
    read = 0
    seen = set()
    safe_makedir("novel")
    with open("novel/hairpin.fa", "w") as fa_handle, open("novel/miRNA.str", "w") as str_handle:
        with open(csv_file) as in_handle:
            for line in in_handle:
                if line.startswith("mature miRBase miRNAs detected by miRDeep2"):
                    break
                if line.startswith("novel miRNAs predicted"):
                    read = 1
                    line = in_handle.next()
                    continue
                if read and line.strip():
                    cols = line.strip().split("\t")
                    name, start, score = cols[0], cols[16], cols[1]
                    if score < 1:
                        continue
                    m5p, m3p, pre = cols[13], cols[14], cols[15].replace('u','t').upper()
                    m5p_start = cols[15].find(m5p) + 1
                    m3p_start = cols[15].find(m3p) + 1
                    m5p_end = m5p_start + len(m5p) - 1
                    m3p_end = m3p_start + len(m3p) - 1
                    if m5p in seen:
                        continue
                    print >>fa_handle, (">{sps}-{name} {start}\n{pre}").format(**locals())
                    print >>str_handle, (">{sps}-{name} ({score}) [{sps}-{name}-5p:{m5p_start}-{m5p_end}] [{sps}-{name}-3p:{m3p_start}-{m3p_end}]").format(**locals())
                    seen.add(m5p)
    return op.abspath("novel")

if __name__ == "__main__":
	out_file = sys.argv[0]
	species = sys.argv[1]
	novel_db = _parse_novel(out_file, species)

