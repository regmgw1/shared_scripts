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
# /peaks_in_segmentation_bedtools.pl
#################################################################

use strict;

unless (@ARGV ==5) {
        die "\n\nUsage:\n ./peaks_in_segmentation_bedtools.pl path2dmrs path2cellList path2segments gff(0)_or_bed(1) path2output\nPlease try again.\n\n\n";}

my $path2dmrs = shift;
my $path2list = shift;
my $path2segments = shift;
my $file = shift;
my $path2output = shift;

my $intersect_thresh = 0.51; #hard-coded as don't want multiple hits - this may need more testing

my (@feats, %mat, %count, %states_count);

my %states = ("CTCF enriched" => "CTCF", "Predicted Enhancer" => "E", "Predicted Promoter Flank" => "PF", "Predicted Promoter with TSS" => "TSS", "Predicted Repressed/Low Activity" => "R", "Predicted Transcribed Region" => "T", "Predicted Weak Enhancer/Cis-reg element" => "WE");

open (OUT, ">$path2output") or die "Can't open $path2output for writing";
print OUT "DMR";
open (COUNT, ">$path2output"."_counts") or die "Can't open $path2output"."_counts for writing";

#for each cell type listed in file $path2list, perform intersectBed vs dmrs. segmentation info for each dmr for each cell type stored in hash %mat
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
				$mat{$coords}{$line} = "NA";
			}
		}
		elsif ($file == 1)
		{
			if ($dmr =~m/chr/)
			{
				$dmr =~s/chr//;
				my @elems = split/\t/, $dmr;
				my $coords = "$elems[0]".":$elems[1]"."-$elems[2]";
				$mat{$coords}{$line} = "NA";
			}
		}
		else
		{
			die "DMR file must be gff (0) or bed (1)!";
		}
		print TMP "$dmr";
		
	}
	close DMR;
	close TMP;
	my $seg_file = "$line"."_segmentation.bed";
	my @count = `intersectBed -a tmp$time.tmp -b $path2segments/$seg_file -wa -wb -f $intersect_thresh`;
	foreach my $out (@count)
	{
		my @elems = split/\t/, $out;
		if ($file == 0)
		{
			my $coords = "$elems[0]".":$elems[3]"."-$elems[4]";
			$mat{$coords}{$line} = $elems[$#elems];
		}
		else
		{
			my $coords = "$elems[0]".":$elems[1]"."-$elems[2]";
			$mat{$coords}{$line} = $elems[$#elems];
		}
		
	}
	my @count2 = `intersectBed -a $path2segments/$seg_file -b tmp$time.tmp -wa -wb -f $intersect_thresh`;
	foreach my $out (@count2)
	{
		my @elems = split/\t/, $out;
		if ($file == 0)
		{
			my $coords = "$elems[4]".":$elems[7]"."-$elems[8]";
			$mat{$coords}{$line} = $elems[3];
		}
		else
		{
			my $coords = "$elems[4]".":$elems[5]"."-$elems[6]";
			$mat{$coords}{$line} = $elems[3];
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
# contents of %mat are walked through and printed out. the segmentation state is converted to abbreviation using %states. Count of each state for each cell type stored in %states_count
foreach my $coord (sort(keys %mat))
{
	print OUT "$coord";
	foreach my $feat (sort(keys %{$mat{$coord}}))
	{
		if ($mat{$coord}{$feat} eq "NA")
		{
			print OUT "\tNA";
		}
		else
		{
			chomp($mat{$coord}{$feat});
			if ($mat{$coord}{$feat} =~m/(.+)-\s$feat/)
			{
				my $tmp = $1;
				chop $tmp;
				#$tmp =~s/\s/_/g;
				print OUT "\t$states{$tmp}";
				if (exists $states_count{$states{$tmp}}{$feat})
				{
					$states_count{$states{$tmp}}{$feat}++;
				}
				else
				{
					$states_count{$states{$tmp}}{$feat} = 1;
				}
			}
			$count{$feat}++;
		}
	}
	print OUT "\n";
}
close OUT; 
# want to output total count of each segment type for each cell type
print COUNT "CellType";
foreach my $seg (sort(keys %states))
{
	print COUNT "\t$states{$seg}";
}
print COUNT "\n";
foreach my $feat (@feats)
{
	print COUNT "$feat";
	foreach my $seg (sort(keys %states))
	{
		if (exists $states_count{$states{$seg}}{$feat})
		{
			print COUNT "\t$states_count{$states{$seg}}{$feat}";
		}
		else
		{
			print COUNT "\t0";
		}
	}
	print COUNT "\n";
}
close COUNT;
#print key to STDOUT
print "\nKEY (more info see http://www.ensembl.org/info/docs/funcgen/regulatory_segmentation.html):\n\n";
foreach my $entry (sort(keys %states))
{
	print "$states{$entry} = $entry\n";
}
