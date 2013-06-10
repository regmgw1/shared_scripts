#!/usr/bin/perl -w

=head2 Authors

=head3 Created by

Gareth Wilson
gareth.wilson@cancer.ucl.ac.uk

=head2 Description
Wrapper script to run peaks_in_features_bedtools.pl as part of the medusa pipeline
=head2 Usage

Usage: ./peaks_in_features_bedtools_wrap.pl path2peakroot path2peaklist path2featureList path2features intersect_thresh path2output

=cut

#################################################################
# peaks_in_features_bedtools_wrap.pl
#################################################################

use strict;

unless (@ARGV ==11) {
        die "\n\nUsage:\n ./peaks_in_features_bedtools_wrap.pl path2peakroot path2peaklist path2featureList path2repeatlist path2features intersect_thresh gff(0)_or_bed(1) path2scripts path2output path2celllist species\nPlease try again.\n\n\n";}

#for pipeline use
my $path2dmrs = shift;
my $path2peaklist = shift;
# for peaks_in_features_bedtools use
my $path2list = shift;
my $path2repeatlist = shift;
my $path2feature = shift;
my $intersect_thresh = shift;
my $dmr_file_type = shift; # 0 for gff, 1 for bed
my $path2scripts = shift;
my $path2output = shift;
my $path2celllist = shift;
my $species = shift;

my (@samples, @dmr_files);

# open list of directory names containing peaks and store in array
open (IN, "$path2peaklist" ) or die "Can't open $path2peaklist for reading";
while (my $line = <IN>)
{
	chomp $line;
	$line =~s/.txt//;
	push @dmr_files, $line;
}
close IN;

foreach my $dmr_file (@dmr_files)
{
	open (OUT, ">$path2output/$dmr_file"."_peak_feature_count.txt") or die "Can't open $path2output/$dmr_file"."_peak_feature_count.txt for writing";
	my $outputPl = "$path2output/$dmr_file"."_peak_feature_matrix.txt";
	my $inputPl = "$path2dmrs/$dmr_file".".bed";
	my @counts = `perl $path2scripts/peaks_in_features_bedtools_v2.pl  $inputPl $path2list $path2feature $intersect_thresh $dmr_file_type $outputPl`;
	foreach my $count (@counts)
	{
		print OUT "$count";
	}
	close OUT;
	open (OUT, ">$path2output/$dmr_file"."_peak_repeat_count.txt") or die "Can't open $path2output/$dmr_file"."_peak_repeat_count.txt for writing";
	$outputPl = "$path2output/$dmr_file"."_peak_repeat_matrix.txt";
	$inputPl = "$path2dmrs/$dmr_file".".bed";
	@counts = `perl $path2scripts/peaks_in_repeats_bedtools.pl  $inputPl $path2repeatlist $path2feature $intersect_thresh $dmr_file_type $outputPl`;
	foreach my $count (@counts)
	{
		print OUT "$count";
	}
	close OUT;
	$species = lc($species);
	if ($species eq "human")
	{
		open (OUT, ">$path2output/segmentation_key.txt") or die "Can't open $path2output/segmentation_key.txt for writing";
		$outputPl = "$path2output/$dmr_file"."_peak_segmentation_matrix.txt";
		$inputPl = "$path2dmrs/$dmr_file".".bed";
		@counts = `perl $path2scripts/peaks_in_segmentation_bedtools.pl  $inputPl $path2celllist $path2feature/regSegments $dmr_file_type $outputPl`;
		foreach my $count (@counts)
		{
			print OUT "$count";
		}
		close OUT;
	}
}
