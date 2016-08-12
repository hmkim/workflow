#!/usr/bin/perl

=head1

 generate small RNA length distribution data and images

 Yi Zheng

=cut

use strict;
use warnings;
use IO::File;

my $usage = qq'
perl $0 input_seq  seq_type (clean, uniq; default = clean)

* clean: the seq ID does not have read number info
* uniq: the seq ID has read number info
* uniq seq ID: ID00001-29 (ID is ID00001, 29 is the count of reads)
';

my $input_seq = shift || die $usage;
my $seq_type = shift || 'clean';

if ($seq_type eq "uniq" || $seq_type eq "clean") {} else { die "Error, seq type is not correct.\n"; }
my $key = $input_seq; $key =~ s/\..*$//;
my $output_table = $key.".table";
my $output_plots = $key.".pdf";
my $output_image = $key.".png";

#################################################################
# Parse the sequence file and get length distribution 		#
#################################################################
my (    $seq_id_info,   # include sequence id and sequence description, separated by spaces
        $seq_id,        # sequence id
        $seq_desc,      # sequence description
        $format,        # format
        $sequence,      # sequence
        $seq_length,    # sequence length
	$max_len,	# maxinum length 
	$min_len,	# mininum length
        $seq_num,       # sequence number
	$uniq_count	# uniq count for uniq reads
);

my %length_dist; $seq_num = 0;

my $fh = IO::File->new($input_seq) || die "Can not open input seq $input_seq $!\n";
while(<$fh>)
{
	chomp;
	$seq_id_info = $_;
	if      ($seq_id_info =~ m/^>/) { $format = 'fasta'; $seq_id_info =~ s/^>//; }
	elsif   ($seq_id_info =~ m/^@/) { $format = 'fastq'; $seq_id_info =~ s/^@//; }
	else    { die "Error at seq format: $seq_id_info\n"; }
	($seq_id, $seq_desc) = split(/\s+/, $seq_id_info, 2);
	unless ($seq_desc) { $seq_desc = ""; }

	$sequence = <$fh>; chomp($sequence);
	$seq_length = length($sequence);

	if ($seq_type eq 'clean') { 
		$seq_num++;
	 } else {
		my @nn = split(/\_/, $seq_id);
		$uniq_count = $nn[scalar(@nn)-1];
        $uniq_count =~ s/x//;
        #print $uniq_count."\n";
		if ($uniq_count < 1) { die "Error, the uniq count for read is not correct $seq_id $uniq_count\n"; }
		$seq_num = $seq_num + $uniq_count;
	}

	if ($max_len) { if ($seq_length > $max_len) { $max_len = $seq_length; }  }
	else { $max_len = $seq_length; }
	if ($min_len) { if ($seq_length < $min_len) { $min_len = $seq_length; }  }
	else { $min_len = $seq_length; }

	if ($seq_type eq 'clean')
	{
		if ( defined $length_dist{$seq_length} ) { $length_dist{$seq_length}++; }
		else { $length_dist{$seq_length} = 1; }
	} 
	else
	{
		if ( defined $length_dist{$seq_length} ) { $length_dist{$seq_length} = $length_dist{$seq_length} + $uniq_count; }
		else { $length_dist{$seq_length} = $uniq_count; }
	}

	if ($format eq 'fastq') { <$fh>; <$fh>; }
}
$fh->close;

my $out = IO::File->new(">".$output_table) || die "Can not open output table file $output_table $!\n";
#print $out "Length\tReads\t%Frequence\n";
foreach my $len (sort keys %length_dist) {
	my $freq = sprintf('%.4f', $length_dist{$len}/$seq_num);
	$freq = $freq * 100;
	print $out "$len\t$length_dist{$len}\t$freq\n";
}
$out->close;

#################################################################
# draw length distribution using R				#
#################################################################

my $R_LD =<< "END";
a<-read.table("$output_table")
x<-a[,1]
y<-a[,2]
dat <- data.frame(fac = rep(x, y))
#pdf("$output_plots",width=12w,height=6)
png("$output_image", width=825, height=600, res=72)
barplot(table(dat)/sum(table(dat)), col="lightblue", xlab="Length(nt)", ylab="Frequency", main="Length distribution")
dev.off()
END

open R,"|/usr/bin/R --vanilla --slave" or die $!;
print R $R_LD;
close R;

#################################################################
# convert pdf file to png format				#
#################################################################		
#my $cmd_convert = "convert $output_plots $output_image";
#system($cmd_convert) && die "Error in command: $cmd_convert\n";

