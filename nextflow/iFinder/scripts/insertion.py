"""
pull out insertion sites from an alignment file aligned with bwa mem
BWA-mem sets the supplementary alignment flag (2048) and sets a SA tag for
reads with supplementary alignments
"""
import os
import pysam
import pybedtools as bedtools
import re
import toolz as tz
from argparse import ArgumentParser

SUPPLEMENTARY_FLAG = 2048
CIGAR_OTHER = 0
CIGAR_SOFT_CLIPPED = 4
CIGAR_HARD_CLIPPED = 5

def convert_cigar_char(char):
    if char == "S":
        return CIGAR_SOFT_CLIPPED
    elif char == "H":
        return CIGAR_HARD_CLIPPED
    else:
        return CIGAR_OTHER

def get_SA_cigar(read):
    SA = [x for x in read.tags if x[0] == 'SA']
    SA = SA[0] if SA else None
    items = SA[1].split(",")
    cigar = items[3]
    bases = [x for x in re.compile("([A-Z])").split(cigar) if x]
    tuples = [(convert_cigar_char(x[1]), x[0]) for x in tz.partition(2, bases)]
    return tuples

def is_supplementary(read):
    return read.flag > SUPPLEMENTARY_FLAG

def get_tag(read, tag):
    tag = [x for x in read.tags if x[0] == tag]
    tag = tag[0] if tag else None
    tag = tag[1] if tag else None
    return tag

def get_SA_items(read):
    SA = [x for x in read.tags if x[0] == 'SA']
    SA = SA[0] if SA else None
    items = SA[1].split(",")
    chrom = items[0]
    pos = int(items[1]) - 1
    strand = items[2]
    cigar = items[3]
    bases = [int(x) for x in re.compile("[A-Z]").split(cigar) if x]
    sa_mapq = items[4].split(";")[0]
    nm = items[5].split(";")[0]
    return chrom, pos, pos + sum(bases), strand, sa_mapq, nm, cigar

def get_SA_tag(read):
    SA = [x for x in read.tags if x[0] == 'SA']
    SA = SA[0] if SA else None
    return SA

def has_supplementary(read):
    return get_SA_tag(read)

def supplementary_contig(read):
    SA = get_SA_tag(read)
    contig = SA[1].split(',')[0]
    return contig

def is_chimera(read, contig, contig_name):
    scontig = supplementary_contig(read)
    if (read.rname == contig) & (scontig != contig_name):
        return True
    elif (read.rname != contig) & (scontig == contig_name):
        return True
    else:
        return False

def parse_cigar_tuples(tuples):
    d = {"clipped_front": 0,
         "insertions": 0,
         "deletions": 0,
         "clipped": 0,
         "matched": 0,
         "other": 0}
    for i, t in enumerate(tuples):
        if t[0] == 4 or t[0] == 5:
            d["clipped"] += t[1]
            if i == 0:
                d['clipped_front'] += t[1]
        elif t[0] == 1:
            d["insertions"] += t[1]
        elif t[0] == 2:
            d["deletions"] += t[1]
        elif t[0] == 0:
            d["matched"] += t[1]
        else:
            d["other"] += t[1]
    return d

def skip_read(read, contig, virus_contig, in_handle, skip_secondary=True):
    skip = False
    if read.is_secondary:
        skip = True
    elif is_supplementary(read) and skip_secondary:
        skip = True
    elif not has_supplementary(read):
        skip = True
    elif not is_chimera(read, contig, virus_contig):
        skip=True
    if read.is_duplicate:
        skip=True
    return skip

def make_hiv_feature_function(bed_file):
    """
    Create a function that takes a chromosome and position and returns all
    overlapping subfeatures from the HIV BED file.

    The BED file must be formatted like this:
    subtype_b       1       634     5-ltr

    calling the returned function with foo("subtype_b", 200) will return 5-ltr
    """
    bed = bedtools.BedTool(bed_file)
    def hiv_overlap(chrom, pos):
        region = bedtools.BedTool("%s %s %s" % (chrom, pos, pos), from_string=True)
        overlap = region.intersect(bed, wao=True)
        features = set([x.fields[6] for x in overlap if x.fields[6] != "."])
        return ",".join(features)
    return hiv_overlap

def make_genomic_feature_function(bed_file):
    """
    Create a function that takes a chromosome and position and returns all
    overlapping subfeatures from the genomic feature BED file.

    The BED file must be formatted like this:
    subtype_b       1       634     5-ltr

    calling the returned function with foo("subtype_b", 200) will return 5-ltr
    """
    bed = bedtools.BedTool(bed_file)
    def overlap_fn(chrom, pos):
        region = bedtools.BedTool("%s %s %s" % (chrom, pos, pos), from_string=True)
        overlap = region.intersect(bed, wao=True)
        symbols = set([x.fields[6] for x in overlap if x.fields[6] != "."])
        symbols = [x for x in symbols if x != "NA"]
        features = set([x.fields[7] for x in overlap if x.fields[7] != "-1"])
        return [",".join(features), ",".join(symbols)]
    return overlap_fn

def igv_reads(bamfile, virus_contig):
    igv_chimeric_file = os.path.splitext(args.bamfile)[0] + ".chimeric.igv.bam"
    if os.path.exists(igv_chimeric_file):
        return igv_chimeric_file
    with pysam.Samfile(args.bamfile, "rb") as in_handle, \
         pysam.Samfile(igv_chimeric_file, "wb", template=in_handle) as out_handle:
        contig = in_handle.gettid(args.virus_contig)
        for read in in_handle:
            if skip_read(read, contig, args.virus_contig, in_handle, False):
                continue
            chrom = in_handle.getrname(read.tid)
            out_handle.write(read)
    return igv_chimeric_file

def keep_chimeric_reads(bamfile, virus_contig):
    chimeric_file = os.path.splitext(args.bamfile)[0] + ".chimeric.bam"
    if os.path.exists(chimeric_file):
        return chimeric_file
    with pysam.Samfile(args.bamfile, "rb") as in_handle, \
         pysam.Samfile(chimeric_file, "wb", template=in_handle) as out_handle:
        contig = in_handle.gettid(args.virus_contig)
        for read in in_handle:
            if skip_read(read, contig, args.virus_contig, in_handle, True):
                continue
            chrom = in_handle.getrname(read.tid)
            out_handle.write(read)
    return chimeric_file

def get_virus_coordinates(read, in_handle, virus_contig):
    chrom = in_handle.getrname(read.tid)
    if chrom == virus_contig:
        return [chrom, read.pos, read.pos + read.reference_length]
    SA_chrom, SA_pos, SA_end, SA_strand, SA_mapq, SA_nm, cigar = get_SA_items(read)
    if SA_chrom == virus_contig:
        return [SA_chrom, SA_pos, SA_end]

def get_human_coordinates(read, in_handle, virus_contig):
    chrom = in_handle.getrname(read.tid)
    if chrom != virus_contig:
        return [chrom, read.pos, read.pos + read.reference_length]
    SA_chrom, SA_pos, SA_end, SA_strand, SA_mapq, SA_nm, cigar = get_SA_items(read)
    if SA_chrom != virus_contig:
        return [SA_chrom, SA_pos, SA_end]

def get_human_cigar(read, in_handle, virus_contig):
    chrom = in_handle.getrname(read.tid)
    if chrom != virus_contig:
        return read.cigar
    SA_chrom, SA_pos, SA_end, SA_strand, SA_mapq, SA_nm, SA_cigar = get_SA_items(read)
    if SA_chrom != virus_contig:
        return get_SA_cigar(read)

def get_clipped_end(read, in_handle, virus_contig):
    cigar = get_human_cigar(read, in_handle, virus_contig)
    first = 0
    last = 0
    if cigar[0][0] == 4 or cigar[0][0] == 5:
        first = cigar[0][1]
    if cigar[-1][0] == 4 or cigar[-1][0] == 5:
        last = cigar[-1][1]
    if first > last:
        return "left"
    elif last > first:
        return "right"
    else:
        return "unknown"

def get_integration_position(read, in_handle, virus_contig):
    clipped_end = get_clipped_end(read, in_handle, virus_contig)
    coords = get_human_coordinates(read, in_handle, virus_contig)
    if clipped_end == "left":
        return coords[1]
    elif clipped_end == "right":
        return coords[2]
    else:
        return None

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--name", default=None, help="Short name of sample.")
    parser.add_argument("--reassign-from", default=None, help="")
    parser.add_argument("--reassign-to", default=None, help="")
    parser.add_argument("--human-bed", help="BED file of human features")
    parser.add_argument("--virus-bed", help="BED file of viral features")
    parser.add_argument("bamfile", help="BAM file to call insertion sites from")
    parser.add_argument("virus_contig", help="Name of virus contig")
    args = parser.parse_args()
    HEADER_FIELDS = ["rid", "first", "chrom", "pos", "end", "strand", "mapq",
                     "SA_chrom", "SA_pos", "SA_end", "SA_strand",
                     "SA_mapq", "SA_nm",
                     "as", "xs", "clipped", "insertions", "deletions",
                     "matched", "clipped_front", "human_chrom", "human_pos",
                     "human_end", "virus_chrom", "virus_pos", "virus_end",
                     "integration"]
    FORMAT = ("{rid}\t{first}\t{chrom}\t{pos}\t{end}\t{strand}\t{mapq}\t{SA_chrom}\t"
              "{SA_pos}\t{SA_end}\t{SA_strand}\t{SA_mapq}\t{SA_nm}\t{AS}\t{XS}\t{clipped}\t"
              "{insertions}\t{deletions}\t{matched}\t{clipped_front}\t{human_chrom}\t"
              "{human_pos}\t{human_end}\t{virus_chrom}\t{virus_pos}\t{virus_end}\t"
              "{integration}")
    OUT_BED = ("{human_chrom}\t{integration}\t{integration_end}\t{rid}\t{name}\t"
               "{virus_pos}\t{orientation}")

    chimeric_file = keep_chimeric_reads(args.bamfile, args.virus_contig)
    igv_file = igv_reads(args.bamfile, args.virus_contig)

    if args.virus_bed:
        FORMAT += "\t{virus_feature}"
        HEADER_FIELDS += ["virus_feature"]
    if args.human_bed:
        FORMAT += "\t{symbol}\t{feature}"
        HEADER_FIELDS += ["symbol", "feature"]
        hiv_feature = make_hiv_feature_function(args.virus_bed)

    if args.human_bed:
        human_feature = make_genomic_feature_function(args.human_bed)

    HEADER = "\t".join(HEADER_FIELDS)
    name = args.name if args.name else os.path.splitext(chimeric_file)[0]

    with pysam.Samfile(chimeric_file, "rb") as in_handle:
        contig = in_handle.gettid(args.virus_contig)
#        print HEADER

        for read in in_handle:
            chrom = in_handle.getrname(read.tid)
            rid = read.qname
            first = read.is_read1
            strand = "-" if read.is_reverse else "+"
            pos = read.pos
            end = read.pos + read.reference_length
            AS = get_tag(read, "AS")
            XS = get_tag(read, "XS")
            mapq = read.mapq
            virus_chrom, virus_pos, virus_end = get_virus_coordinates(read, in_handle, args.virus_contig)
            human_chrom, human_pos, human_end = get_human_coordinates(read, in_handle, args.virus_contig)
            SA_chrom, SA_pos, SA_end, SA_strand, SA_mapq, SA_nm, SA_cigar = get_SA_items(read)
            cigar = parse_cigar_tuples(read.cigar)
            clipped = cigar["clipped"]
            insertions = cigar["insertions"]
            deletions = cigar["deletions"]
            matched = cigar["matched"]
            clipped_front = cigar["clipped_front"]
            integration = get_integration_position(read, in_handle, args.virus_contig)
            if not integration:
                continue
	        integration_end = integration + 1
            orientation = get_clipped_end(read, in_handle, args.virus_contig)
            if args.virus_bed:
                virus_feature = hiv_feature(*get_virus_coordinates(read, in_handle, args.virus_contig))
            if args.human_bed:
                symbol, feature = human_feature(*get_human_coordinates(read, in_handle, args.virus_contig))
            print OUT_BED.format(**locals())
#            print FORMAT.format(**locals())
