#!/usr/bin/perl -w

=head2 Authors

=head3 Created by

Gareth Wilson
gareth.wilson@cancer.ucl.ac.uk

=head2 Description
Use bedTools 'intersectBed' to determine how many reads are found under peaks. Quicker than peak_explorer.pl.
=head2 Usage

Usage: ./peaks_explorer_bedtools.pl sample_list path2bedDir path2peakroot path2peaklist path2output

=cut

#################################################################
# peaks_explorer.pl
#################################################################

use strict;

unless (@ARGV ==5) {
        die "\n\nUsage:\n ./peaks_explorer_bedtools.pl sample_list path2bedDir path2peakroot path2peaklist path2output\nPlease try again.\n\n\n";}

my $path2samplelist = shift;
my $path2bed = shift;
my $path2peakroot = shift;
my $path2peaklist = shift;
my $path2output = shift;

my (@samples, @peak_subs);

# open list of sample names and store in array
open (IN, "$path2samplelist" ) or die "Can't open $path2samplelist for reading";
while (my $line = <IN>)
{
	chomp $line;
	push @samples, $line;
}
close IN;

# open list of directory names containing peaks and store in array
open (IN, "$path2peaklist" ) or die "Can't open $path2peaklist for reading";
while (my $line = <IN>)
{
	chomp $line;
	chop $line;
	push @peak_subs, $line;
}
close IN;

foreach my $peaksub (@peak_subs)
{
	open (OUT, ">$path2output/$peaksub"."_peak_explore.txt") or die "Can't open $path2output/$peaksub"."_peak_explore.txt for writing";
	my $peakfile = $peaksub;
	my (%peak, %peak_pos, %peak_neg, %coord);
	$peakfile =~s/.*Binary/binary/;
	print OUT "Chr\tPeakStart\tPeakStop";
	foreach my $sample (sort @samples)
	{
		print "$sample\n";
		print OUT "\t$sample\t$sample Pos\t$sample Neg";
		my @overlaps = `intersectBed -a $path2peakroot/$peaksub/$peakfile"".gff -b $path2bed/$sample"".bed -wao`;
		foreach my $overlap (@overlaps)
		{
			chomp $overlap;
			my @elems = split /\t/, $overlap;
			$coord{$elems[3]} = $elems[4];
			if (exists $peak{$elems[0]}{$elems[3]}{$sample})
			{
				$peak{$elems[0]}{$elems[3]}{$sample}++;
				if ($elems[14] eq "+")
				{
					$peak_pos{$elems[0]}{$elems[3]}{$sample}++;
				}
				else
				{
					$peak_neg{$elems[0]}{$elems[3]}{$sample}++;
				}
			}
			else
			{
				if ($elems[15] == 0)
				{
					$peak{$elems[0]}{$elems[3]}{$sample} = 0;
				}
				else
				{
					$peak{$elems[0]}{$elems[3]}{$sample} = 1;
					if ($elems[14] eq "+")
					{
						$peak_pos{$elems[0]}{$elems[3]}{$sample} = 1;
					}
					else
					{
						$peak_neg{$elems[0]}{$elems[3]}{$sample} = 1;
					}
				}
			}
		}
		
	}
	print OUT "\n";
	foreach my $chr_out (sort keys %peak)
	{
		for my $begin_out (sort numerically keys %{$peak{$chr_out}})
		{
			print OUT "$chr_out\t$begin_out\t$coord{$begin_out}";
			foreach my $samp_out (sort keys %{$peak{$chr_out}{$begin_out}})
			{
				print OUT "\t$peak{$chr_out}{$begin_out}{$samp_out}";
				if (exists $peak_pos{$chr_out}{$begin_out}{$samp_out})
				{
					print OUT "\t$peak_pos{$chr_out}{$begin_out}{$samp_out}";
				}
				else
				{
					print OUT "\t0";
				}
				if (exists $peak_neg{$chr_out}{$begin_out}{$samp_out})
				{
					print OUT "\t$peak_neg{$chr_out}{$begin_out}{$samp_out}";
				}
				else
				{
					print OUT "\t0";
				}								
			}
			print OUT "\n";
		}
	}
	close OUT;
}
sub numerically {$a<=>$b};				
