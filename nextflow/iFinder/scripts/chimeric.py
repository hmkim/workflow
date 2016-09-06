"""
find all chimeric reads chimeric for the virus specified as virus_contig
from a BAM file mapped with bwa-mem with the -Y flag set, to stop hard
clipping the chimeric alignments
reports both the primary and secondary alignment
"""

import pysam
import os
from argparse import ArgumentParser
SUPPLEMENTARY_FLAG = 2048

def get_SA_tag(read):
    SA = [x for x in read.tags if x[0] == 'SA']
    SA = SA[0] if SA else None
    return SA

def supplementary_contig(read):
    SA = get_SA_tag(read)
    if not SA:
        return None
    contig = SA[1].split(',')[0]
    return contig

def is_chimera(read, chrom, virus_contig):
    schrom = supplementary_contig(read)
    if (chrom == virus_contig) & (schrom != virus_contig):
        return True
    elif (chrom != virus_contig) & (schrom == virus_contig):
        return True
    else:
        return False

def is_supplementary(read):
    return read.flag > SUPPLEMENTARY_FLAG

def is_chimeric_read(read, chrom, virus_contig):
    is_supplementary_virus = is_supplementary(read) and chrom == virus_contig
    return is_supplementary_virus or is_chimera(read, chrom, virus_contig)

def chimeric_reads(bamfile, virus_contig, duplicates):
    igv_chimeric_file = os.path.splitext(args.bamfile)[0] + ".chimeric.igv.bam"
    if os.path.exists(igv_chimeric_file):
        return igv_chimeric_file
    with pysam.Samfile(args.bamfile, "rb") as in_handle, \
         pysam.Samfile(igv_chimeric_file, "wb", template=in_handle) as out_handle:
        for read in in_handle:
            try:
                chrom = in_handle.getrname(read.tid)
            except ValueError:
                continue
            if not is_chimeric_read(read, chrom, virus_contig):
                continue
            if read.qname in duplicates:
                continue
            out_handle.write(read)
    return igv_chimeric_file

def get_duplicates(bam_file):
    s = set()
    with pysam.Samfile(bam_file, "rb") as in_handle:
        for read in in_handle:
            if read.is_duplicate:
                s.update([read.qname])
    return s

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("virus_contig", help="Name of virus contig")
    parser.add_argument("bamfile", help="BAM file")
    args = parser.parse_args()
    duplicates = get_duplicates(args.bamfile)
    chimeric_reads(args.bamfile, args.virus_contig, duplicates)
