"""
given a BAM file of chimeric alignments both primary and secondary, call

a) the insertion site on the human genome
b) the orienation of the HIV sequence (attached to 3' or 5' end)
"""

from argparse import ArgumentParser
import pysam
import toolz as tz
import re

CIGAR_SOFT_CLIP = 4
CIGAR_HARD_CLIP = 5

def get_virus_end(read):
    cigar = read.cigartuples
    first = 0
    last = 0
    if not cigar:
        return "unknown"
    if cigar[0][0] == CIGAR_SOFT_CLIP or cigar[0][0] == CIGAR_HARD_CLIP:
        first = int(cigar[0][1])
    if cigar[-1][0] == CIGAR_SOFT_CLIP or cigar[-1][0] == CIGAR_HARD_CLIP:
        last = int(cigar[-1][1])
    # if both ends are clipped > 20 bases, skip this alignment
    if (first > 20) and (last > 20):
        return "unknown"
    if first > last:
        return "R"
    elif last > first:
        return "L"
    else:
        return "unknown"

def get_insertion_end(cigar):
    first = 0
    last = 0
    if cigar[0][0] == CIGAR_SOFT_CLIP or cigar[0][0] == CIGAR_HARD_CLIP:
        first = int(cigar[0][1])
    if cigar[-1][0] == CIGAR_SOFT_CLIP or cigar[-1][0] == CIGAR_HARD_CLIP:
        last = int(cigar[-1][1])
    # if both ends are clipped > 20 bases, skip this alignment
    if (first > 20) and (last > 20):
        return "unknown"
    if first > last:
        return "L"
    elif last > first:
        return "R"
    else:
        return "unknown"

def get_virus_code(read):
    read_direction = "For" if read.is_read1 else "Rev"
    virus_end = get_virus_end(read)
    strand = "-" if read.is_reverse else "+"
    return read_direction + virus_end + strand

def get_integration_code(read, insertion_end, strand):
    read_direction = "For" if read.is_read1 else "Rev"
    return read_direction + insertion_end + strand

def call_orientation(virus_code):
    if virus_code in ["ForR+", "RevR-", "RevR+", "ForR-"]:
        return "5prime"
    elif virus_code in ["RevL-", "ForL-", "RevL+", "ForL+"]:
        return "3prime"
    else:
        return "unknown"

def total_bases_cigar(cigar):
    bases = [x for x in re.compile("([A-Z])").split(cigar) if x]
    z = [count_tuple(x) for x in tz.partition(2, bases)]
    return z

def count_tuple(tp):
    if tp[1] in ["S", "H"]:
        return 0
    elif tp[1] in ["D"]:
        return -int(tp[0])
    else:
        return int(tp[0])

def get_SA_items(read):
    SA = [x for x in read.tags if x[0] == 'SA']
    SA = SA[0] if SA else None
    if not SA:
        return None, None, None, None, None, None, None
    items = SA[1].split(",")
    chrom = items[0]
    pos = int(items[1]) - 1
    strand = items[2]
    cigar = items[3]
    total_bases = total_bases_cigar(cigar)
    sa_mapq = items[4].split(";")[0]
    nm = items[5].split(";")[0]
    return chrom, pos, pos + sum(total_bases), strand, sa_mapq, nm, cigar

def get_SA_cigar(read):
    SA = [x for x in read.tags if x[0] == 'SA']
    SA = SA[0] if SA else None
    items = SA[1].split(",")
    cigar = items[3]
    bases = [x for x in re.compile("([A-Z])").split(cigar) if x]
    tuples = [(convert_cigar_char(x[1]), x[0]) for x in tz.partition(2, bases)]
    return tuples

def convert_cigar_char(char):
    CIGAR_OTHER = 0
    CIGAR_SOFT_CLIPPED = 4
    CIGAR_HARD_CLIPPED = 5
    if char == "S":
        return CIGAR_SOFT_CLIPPED
    elif char == "H":
        return CIGAR_HARD_CLIPPED
    else:
        return CIGAR_OTHER

def get_duplicates(bam_file):
    s = set()
    with pysam.Samfile(bam_file, "rb") as in_handle:
        for read in in_handle:
            if read.is_duplicate:
                s.update([read.qname])
    return s

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("bamfile", help="BAM file of chimeric reads")
    parser.add_argument("virus_contig", help="Name of contig of virus.")
    args = parser.parse_args()

    print "read_name chrom pos strand orientation insertion_end mapq virus_pos virus_end seq code seqcode"

    duplicates = get_duplicates(args.bamfile)

    with pysam.Samfile(args.bamfile, "rb") as in_handle:
        count = 0
        for read in in_handle:
            chrom = in_handle.getrname(read.tid)
            if read.is_duplicate:
                continue
            if chrom != args.virus_contig:
                continue
            if read.qname in duplicates:
                continue
            code = get_virus_code(read)
            orientation =  call_orientation(code)
            SA_chrom, SA_pos, SA_end, SA_strand, SA_mapq, SA_nm, SA_cigar = get_SA_items(read)
            if SA_chrom == args.virus_contig:
                continue
            if not SA_chrom:
                continue
            SA_cigar_tuples = get_SA_cigar(read)
            insertion_end = get_insertion_end(SA_cigar_tuples)
            virus_end = get_virus_end(read)
            if insertion_end == "unknown" or virus_end == "unknown":
                continue
            insertion_pos = SA_pos if insertion_end == "L" else SA_end
            end_pos = SA_pos if insertion_end == "R" else SA_end
            icode = get_integration_code(read, insertion_end, SA_strand)
            iorientation = call_orientation(icode)
            if virus_end == "L":
                seqcode = "{orientation}(HIV)-{insertion_pos}-{iorientation}(human)".format(**locals())
            else:
                seqcode = "{iorientation}(human)-{insertion_pos}-{orientation}(HIV)".format(**locals())
            print read.qname, SA_chrom, insertion_pos, SA_strand, orientation, insertion_end, SA_mapq, read.pos, read.aend, read.seq, code, seqcode
