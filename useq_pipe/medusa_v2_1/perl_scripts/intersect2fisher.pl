#!/usr/bin/perl -w

=head2 Authors

=head3 Created by

Gareth Wilson
gareth.wilson@cancer.ucl.ac.uk

#adapted by Anna Koeferle
#input files need to consist of  3 columns (chr, start, end)


=head2 Description
Use bedTools to generate file for input into R for chi-square analysis
=head2 Usage

Usage: ./intersect2fisher.pl path2allPossiblePeaks path2dmrs path2features path2featureList

=cut

#################################################################
# intersect2chi.pl
#################################################################

use strict;

unless (@ARGV == 7) {
        die "\n\nUsage:\n ./intersect2fisher.pl path2dmrs path2dmrFile path2allPossiblePeaks path2features path2featureList path2Rscripts path2output\nPlease try again.\n\n\n";}

my $path2dmrs = shift;
my $path2dmrlist = shift;
my $path2all = shift;
my $path2features = shift;
my $path2list =shift;
my $path2Rscripts = shift;
my $path2output = shift;

my (@dmr_files);
open (IN, "$path2dmrlist" ) or die "Can't open $path2dmrlist for reading";
while (my $line = <IN>)
{
	chomp $line;
	$line =~s/.txt//;
	push @dmr_files, $line;
}
close IN;

foreach my $dmr_set (@dmr_files)
{
	my $time = time();
	open (TMP, ">dmr$time.tmp") or die "Can't open tmp.tmp for writing";
	open (DMR, "$path2dmrs/$dmr_set.bed" ) or die "Can't open $path2dmrs/$dmr_set.bed for reading";
	while (my $dmrt = <DMR>)
	{
		$dmrt =~s/chr//;
		print TMP "$dmrt";
	}
	close DMR;
	close TMP;

	open (TMP, ">peaks$time.tmp") or die "Can't open tmp.tmp for writing";
	open (ALL, $path2all ) or die "Can't open $path2all for reading";
	while (my $peak = <ALL>)
	{       
	        $peak =~s/chr//;
	        print TMP "$peak";
	}
	close ALL;
	close TMP;
	open (OUT, ">$path2output/$dmr_set"."_enrichment.txt") or die "Can't open $path2output/$dmr_set"."_enrichment.txt for writing";
	print OUT "Features\tDMRs_in_feature\tTotal_DMRs\tPotential_in_feature\tTotal_potential\n";
	open (IN, "$path2list" ) or die "Can't open $path2list for reading";
	while (my $line = <IN>)
	{
		chomp $line;
		my $dmr = `bash -c 'cat <(intersectBed -f 0.51 -a dmr$time.tmp -b $path2features/$line/$line.gff -wa | cut -f 1-3 | intersectBed -f 0.51 -a stdin -b peaks$time.tmp -wa -wb | cut -f 4-6) <(intersectBed -f 0.51 -b dmr$time.tmp -a $path2features/$line/$line.gff -wa -wb | cut -f 10-12 |intersectBed -f 0.51 -a stdin -b peaks$time.tmp -wa -wb | cut -f 4-6) <(intersectBed -f 0.51 -a dmr$time.tmp -b $path2features/$line/$line.gff -wa | cut -f 1-3 | intersectBed -f 0.51 -b stdin -a peaks$time.tmp -wa | cut -f 1-3) <(intersectBed -f 0.51 -b dmr$time.tmp -a $path2features/$line/$line.gff -wa -wb | cut -f 10-12 |intersectBed -f 0.51 -b stdin -a peaks$time.tmp -wa | cut -f 1-3) | sortBed | mergeBed | wc -l'`;
		chomp $dmr;
		my $totalfeat = `bash -c 'cat <(intersectBed -f 0.51 -a peaks$time.tmp -b $path2features/$line/$line.gff  -wa | cut -f 1-3) <(intersectBed -f 0.51 -b peaks$time.tmp -a $path2features/$line/$line.gff -wa -wb | cut -f 10-12) | sortBed | mergeBed |  wc -l'`;
		chomp $totalfeat;
		my $totaldmr = `bash -c 'cat <(intersectBed -f 0.51  -a peaks$time.tmp -b dmr$time.tmp -wa | cut -f 1-3) <(intersectBed -f 0.51 -b peaks$time.tmp -a dmr$time.tmp -wa -wb | cut -f 5-7) | sortBed | mergeBed | wc -l'`;
		chomp $totaldmr;
		my $total = `grep ^chr $path2all|wc -l`;
		chomp $total;
		print OUT "$line\t$dmr\t$totaldmr\t$totalfeat\t$total\n";
	}
	close IN;
	close OUT;
	unlink("peaks$time.tmp","dmr$time.tmp");
	system "/$path2Rscripts/features_fisher.R $path2output/ $dmr_set";
}
