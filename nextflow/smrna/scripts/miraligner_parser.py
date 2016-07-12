"3""Functions to access json with classes"""
from collections import defaultdict
import pandas as pd
import sys
import os.path as op


try:
    from seqcluster import prepare_data as prepare
    from seqcluster import templates as template_seqcluster
    from seqcluster.seqbuster import _create_counts, _read_miraligner, _tab_output
except ImportError:
    pass

#sys.path.apeend('/BiO/BioTools/bcbio/data/anaconda/lib/python2.7/site-packages/bcbio')

from bcbio.utils import file_exists, safe_makedir, move_safe, append_stem


class realign:

    def __init__(self):
        self.sequence = ""
        self.precursors = defaultdict(isomir)
        self.score = []
        self.best_hits = [] # maybe sam object?

    def set_precursor(self, precursor, isomir):
        self.precursors[precursor] = isomir

    def remove_precursor(self, precursor):
        del self.precursors[precursor]

class isomir:

    def __init__(self):
        self.t5 = "0"
        self.t3 = "0"
        self.add = []
        self.subs = []
        self.align = None
        self.start = 0
        self.mirna = None

    def format(self, sep="\t"):
        subs = "".join(["".join(map(str, mism)) for mism in self.subs])
        if not subs:
            subs = "0"
        add = "0" if not self.add else self.add
        return "%s%s%s%s%s%s%s" % (subs, sep, add, sep,
                                   self.t5, sep, self.t3)

    def get_score(self, sc):
        for a in self.add:
            if a in ['A', 'T']:
                sc -= 0.25
            else:
                sc -= 0.75
        for e in self.subs:
            sc -= 1
        return sc

def _get_freq(name):
    """
    Check if name read contains counts (_xNumber)
    """
    try:
        counts = int(name.split("_x")[1])
    except:
        return 0
    return counts


def _tab_output(reads, out_file, sample):
    seen = set()
    lines = []
    lines_pre = []
    seen_ann = {}
    dt = None
    with open(out_file, 'w') as out_handle:
        print >>out_handle, "name\tseq\tfreq\tchrom\tstart\tend\tsubs\tadd\tt5\tt3\ts5\ts3\tDB\tprecursor\thits"
        for r, read in reads.iteritems():
            hits = set()
            [hits.add(mature.mirna) for mature in read.precursors.values() if mature.mirna]
            hits = len(hits)
            for p, iso in read.precursors.iteritems():
                if len(iso.subs) > 3 or not iso.mirna:
                    continue
                if (r, iso.mirna) not in seen:
                    seen.add((r, iso.mirna))
                    chrom = iso.mirna
                    if not chrom:
                        chrom = p
                    count = _get_freq(r)
                    seq = reads[r].sequence
                    if iso.get_score(len(seq)) < 1:
                        continue
                    if iso.subs:
                        iso.subs = [] if "N" in iso.subs[0] else iso.subs
                    annotation = "%s:%s" % (chrom, iso.format(":"))
                    res = ("{seq}\t{r}\t{count}\t{chrom}\tNA\tNA\t{format}\tNA\tNA\tmiRNA\t{p}\t{hits}").format(format=iso.format().replace("NA", "0"), **locals())
                    if annotation in seen_ann and seq.find("N") < 0 and seen_ann[annotation].split("\t")[0].find("N") < 0:
                        raise ValueError("Same isomir %s from different sequence: \n%s and \n%s" % (annotation, res, seen_ann[annotation]))
                    seen_ann[annotation] = res
                    lines.append([annotation, chrom, count, sample, hits])
                    lines_pre.append([annotation, chrom, p, count, sample, hits])
                    print >>out_handle, res

    if lines:
        dt = pd.DataFrame(lines)
        dt.columns = ["isomir", "chrom", "counts", "sample", "hits"]
        dt = dt[dt['hits']>0]
        dt = dt.loc[:, "isomir":"sample"]
        dt = dt.groupby(['isomir', 'chrom', 'sample'], as_index=False).sum()
        dt.to_csv(out_file + "_summary")
        dt_pre = pd.DataFrame(lines_pre)
        dt_pre.columns = ["isomir", "mature", "chrom", "counts", "sample", "hits"]
        dt_pre = dt_pre[dt_pre['hits']==1]
        dt_pre = dt_pre.loc[:, "isomir":"sample"]
        dt_pre = dt_pre.groupby(['isomir', 'chrom', 'mature', 'sample'], as_index=False).sum()
        return out_file, dt, dt_pre
    return None

def _read_miraligner(fn):
    """Read ouput of miraligner and create compatible output."""
    reads = defaultdict(realign)
    with open(fn) as in_handle:
        in_handle.next()
        for line in in_handle:
            cols = line.strip().split("\t")
            iso = isomir()
            query_name, seq = cols[1], cols[0]
            chrom, reference_start = cols[-2], cols[3]
            iso.mirna = cols[3]
            subs, add, iso.t5, iso.t3 = cols[6:10]
            if query_name not in reads:
                reads[query_name].sequence = seq
            iso.align = line
            iso.start = reference_start
            iso.subs, iso.add = _parse_mut(subs), add
            #logger.debug("%s %s %s %s %s" % (query_name, reference_start, chrom, iso.subs, iso.add))
            reads[query_name].set_precursor(chrom, iso)
    return reads

def _parse_mut(subs):
    """
    Parse mutation tag from miraligner output
    """
    if subs!="0":
        subs = [[subs.replace(subs[-2:], ""),subs[-2], subs[-1]]]
    return subs

if __name__ == "__main__":
	infile = sys.argv[1]
	out_file = sys.argv[2]
	sample = sys.argv[3]

	reads = _read_miraligner(infile)

	out_dts = []

	out_file, dt, dt_pre= _tab_output(reads, out_file, sample)
	out_dts.append(dt)

	out_dir = "./"
	stem=""
	
	out_novel_isomir = append_stem(op.join(out_dir, "counts.tsv"), stem)
	out_novel_mirna = append_stem(op.join(out_dir, "counts_mirna.tsv"), stem)

	out_files = _create_counts(out_dts, out_dir)
        out_files = [move_safe(out_files[0], out_novel_isomir), move_safe(out_files[1], out_novel_mirna)]

	print out_files

