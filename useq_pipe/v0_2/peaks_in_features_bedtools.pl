#!/usr/bin/perl -w

=head2 Authors

=head3 Created by

Gareth Wilson
gareth.wilson@cancer.ucl.ac.uk

=head2 Description
Utilises bedTools to find overlaps between feature types. Dmr file can be bed or gff. Feature file must be gff. Matrix output is binary i.e. 0 no coverage, 1 dmr covers feature/features.
=head2 Usage

Usage: /peaks_in_features_bedtools.pl path2dmrs path2featureList path2features intersect_threshold path2output

=cut

#################################################################
# /peaks_in_features_bedtools.pl
#################################################################

use strict;

unless (@ARGV ==6) {
        die "\n\nUsage:\n ./peaks_in_features_bedtools.pl path2dmrs path2featureList path2features intersect_threshold gff(0)_or_bed(1) path2output\nPlease try again.\n\n\n";}

my $path2dmrs = shift;
my $path2list = shift;
my $path2feature = shift;
my $intersect_thresh = shift;
my $file = shift;
my $path2output = shift;

my (@feats, %mat, %count);

open (OUT, ">$path2output") or die "Can't open $path2output for writing";
print OUT "DMR";

open (IN, "$path2list" ) or die "Can't open $path2list for reading";
while (my $line = <IN>)
{
	print STDERR $line;
	chomp $line;
	my $time = time();
	open (TMP, ">tmp$time.tmp") or die "Can't open tmp.tmp for writing";
	open (DMR, $path2dmrs ) or die "Can't open $path2dmrs for reading";
	while (my $dmr = <DMR>)
	{
		if ($file == 0)
		{
			if ($dmr =~m/chr/)
			{
				$dmr =~s/chr//;
				my @elems = split/\t/, $dmr;
				my $coords = "$elems[0]".":$elems[3]"."-$elems[4]";
				$mat{$coords}{$line} = 0;
			}
		}
		elsif ($file == 1)
		{
			my @elems = split/\t/, $dmr;
			my $coords = "$elems[0]".":$elems[1]"."-$elems[2]";
			$mat{$coords}{$line} = 0;
		}
		else
		{
			die "DMR file must be gff (0) or bed (1)!";
		}
		print TMP "$dmr";
		
	}
	close DMR;
	close TMP;
	my @count = `intersectBed -a tmp$time.tmp -b $path2feature/$line/$line"."gff -wa -wb -f $intersect_thresh`;
	foreach my $out (@count)
	{
		my @elems = split/\t/, $out;
		if ($file == 0)
		{
			my $coords = "$elems[0]".":$elems[3]"."-$elems[4]";
			$mat{$coords}{$line}++;
		}
		else
		{
			my $coords = "$elems[0]".":$elems[1]"."-$elems[2]";
			$mat{$coords}{$line}++;
		}
		
	}
	my @count2 = `intersectBed -a $path2feature/$line/$line"."gff -b tmp$time.tmp -wa -wb -f $intersect_thresh`;
	foreach my $out (@count2)
	{
		my @elems = split/\t/, $out;
		if ($file == 0)
		{
			my $coords = "$elems[9]".":$elems[12]"."-$elems[13]";
			$mat{$coords}{$line}++;
		}
		else
		{
			my $coords = "$elems[9]".":$elems[10]"."-$elems[11]";
			$mat{$coords}{$line}++;
		}
		
	}
	unlink ("tmp$time.tmp");
	undef(@count);
	push @feats, $line;
}
close IN;

@feats = sort(@feats);
foreach my $feat (@feats)
{
	print OUT "\t$feat";
	$count{$feat} = 0;
}
print OUT "\n";
foreach my $coord (sort(keys %mat))
{
	print OUT "$coord";
	foreach my $feat (sort(keys %{$mat{$coord}}))
	{
		if ($mat{$coord}{$feat} > 0)
		{
			print OUT "\t1";
			$count{$feat}++;
		}
		else
		{
			print OUT "\t0";
		}

	}
	print OUT "\n";
}
close OUT; 
foreach my $feat (@feats)
{
	print "$feat\t$count{$feat}\n";
}
