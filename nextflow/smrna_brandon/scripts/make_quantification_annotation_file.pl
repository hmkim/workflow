use warnings;
use strict;

use Data::Dumper;

if (@ARGV != 2){
	print "Usage: perl $0 <gene_count_file> <in_gtf_file>\n";
	exit;
}

my $gene_count_file = $ARGV[0];
#my $gene_count_file = "/BiO/BioProjects/Nexbio-SmallRNAseq-2016-01/Result/tmp/HTseq_count/TN1509R1166.gencode.genecount.txt";
checkFile($gene_count_file);

my %gene_counts;
read_gene_count_file($gene_count_file, \%gene_counts);

#my $in_gtf_file = "/BiO/BioResources/ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.chrRemove.gtf";
#my $in_gtf_file = "/BiO/BioResources/References/Bos_taurus/UMD3.1.83/Bos_taurus.UMD3.1.83.gtf";
my $in_gtf_file = $ARGV[1];
checkFile($in_gtf_file);

read_gtf_file($in_gtf_file);

sub read_gene_count_file{
	my ($file, $hash_ref) = @_;

	open my $fh, '<:encoding(UTF-8)', $file or die;
	while (my $row = <$fh>) {
		chomp $row;
		my ($gene,$count) = split /\t/, $row;
		$hash_ref->{$gene} = $count;
	}
	close($fh);
	
}

sub read_gtf_file{
	my $file = shift;

	#open my $fh, '<:encoding(UTF-8)', $file or die;
	my $fh;
	if ($file =~ /.gz$/) {
		open $fh, "gunzip -c $file |" || die "can’t open pipe to $file";
	}
	else {
		open $fh, $file || die "can’t open $file";
	}

	while (my $row = <$fh>) {
		chomp $row;
		if ($row =~ /^#/){
			next;
		}
		my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split /\t/, $row;
		my @arr_attribute = split /\;/, $attribute;


		my %attr = attr_to_hash(\@arr_attribute);
		if ($feature eq "gene"){
			my $gene_id = $attr{gene_id};

			my $gene_type = $attr{gene_type};
			my $gene_biotype = $attr{gene_biotype};

			my $gene_source = $attr{gene_source};
			my $gene_status = $attr{gene_status};
			my $gene_name = $attr{gene_name};

			if (!$gene_id and !$gene_name){
				die "ERROR! not defined column <$row>\n";
			}else{
				if (!$gene_name){
					$gene_name = $gene_id;
				}
			}
			if (!$gene_status){
				if ($gene_source){ # for ensembl gtf
					$gene_status = $gene_source
				}else{
					die "ERROR! not defined gene_status <$row>\n";
				}
			}
			if ($gene_biotype){
				$gene_type = $gene_biotype;
			}
			if (!$gene_type and !$gene_biotype){
				die "ERROR! not defined type (gene_type and gene_biotype is no exist) <$row>\n";
			}

			if (defined $gene_counts{$gene_id}){
				my $gene_count = $gene_counts{$gene_id};
				print "$gene_id\t$gene_count\t$gene_id\t$gene_type\t$gene_status\t$gene_name\n";
			}	
		}	
	}
	close($fh);

}
#cat /BiO/BioResources/ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.chrRemove.gtf | awk -F"\t" '($3 == "gene"){print}' | awk -F "\"" '{print $2"\t"$6"\t"$8"\t"$10}' | sort -k 1 | paste /BiO/BioProjects/Nexbio-SmallRNAseq-2016-01/Result/tmp/HTseq_count/TN1509R1166.gencode.genecount.txt - | awk '($1 == $3){print}' > /BiO/BioProjects/Nexbio-SmallRNAseq-2016-01/Result/tmp/HTseq_count/TN1509R1166.gencode.genecount.annotation.txt

sub checkFile{
	my $file = shift;
	if (!-f $file){
		die "ERROR! not found <$file>\n";
	}
}

sub attr_to_hash{
	my $arr_ref = shift;

	my %hash;
	foreach (@$arr_ref){
		$_ = trim($_);
		my ($key, $value) = split / /, $_;
		if ($value =~ /\"/){ $value =~ s/\"//g; }
		$hash{$key} = $value;
	}

	return %hash;	
}

sub trim {
	my @result = @_;

	foreach (@result) {
		s/^\s+//;
		s/\s+$//;
	}

	return wantarray ? @result : $result[0];
}
#seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
#source - name of the program that generated this feature, or the data source (database or project name)
#feature - feature type name, e.g. Gene, Variation, Similarity
#start - Start position of the feature, with sequence numbering starting at 1.
#end - End position of the feature, with sequence numbering starting at 1.
#score - A floating point value.
#strand - defined as + (forward) or - (reverse).
#frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
#attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.
#
