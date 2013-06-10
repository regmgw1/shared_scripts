#!/usr/bin/perl -w

#!/usr/bin/perl -w

=head2 Authors

=head3 Created by

Gareth Wilson
gareth.wilson@cancer.ucl.ac.uk

=head2 Description
to be run on filtered sam files. convert to bed/gff with only unique reads included
=head2 Usage

Usage: ./sam2bedgff.pl path2sam sample species read_length single_end? bed? path2output

=cut

#################################################################
# sam2bedgff.pl - to be run on filtered sam files. convert to bed/gff with only unique reads included
#################################################################

use strict;
use File::Basename;
use IO::File;

unless (@ARGV ==8 ) {
        die "\n\nUsage:\n ./sam2bedgff.pl path2sam sample species read_length chrom_list single_end? bed? path2output\nPlease try again.\n\n\n";}

my $path2sam = shift;
my $sample = shift;
my $species = shift;
my $read_length = shift;
my $chrom_list = shift;
my $se = shift;
my $bed = shift;
my $path2output = shift;

#Determining number of chromosomes based on species
my (@chroms,%typeH);

if ($species eq "single")
{
	@chroms = (1);
}
else
{
	open (IN, "$chrom_list" ) or die "Can't open $chrom_list for reading";
	while (my $line = <IN>)
	{
		chomp $line;
		my @elems = split/\t/, $line;
		my $chrCol = $elems[0];
		if ($chrCol !~m/_/)
		{
			if ($chrCol !~m/^chr/)
			{
				$chrCol = "chr".$chrCol;
			}
			#$chrCol =~s/chr//;
			push @chroms, $chrCol;
			my $fh = new IO::File(">$path2output"."_$chrCol".".sam") or die "Can't open $path2output"."_$chrCol".".sam";
			$typeH{$chrCol} = $fh;
			print "IN $chrCol\n";
		}
		else
		{
			print "EXCLUDED $chrCol - not a major chromosome and will be removed from further analysis\n";
		}
	}
	close IN;
}

my $total_rep = 0;
my $total_read = 0;
my $failed_qc = 0;
my $count_lines = 0;

open (OUT ,">$path2output") or die "Can't open $path2output for writing";

#Separate chromosomes into individual temporary files
#my $chrTemp;
#$chrTemp = $path2output;
#$chrTemp =~s/\.(bed|gff)/\_chr/;

print "Creating temporary chromosome files\n";
if ($species ne "single")
{
	open (IN, "$path2sam" ) or die "Can't open $path2sam for reading";
	while (my $line = <IN>)
	{
		my @elems = split/\t/,$line;
		my $chrom = $elems[2];
		if ($chrom !~m/_/)
		{
			if ($chrom !~m/^chr/)
			{
				$chrom = "chr".$chrom;
			}
			#$chrom =~s/^chr//;
			my $fh = $typeH{$chrom};
			print $fh "$line";
		}
	}
	close IN;
}
# close all output files
foreach my $fh (values (%typeH))
{
	close $fh;
}
print "Creating .bed file\n";
foreach my $chr (@chroms)
{
	my $count_rep = 0;
	my $count_reads = 0;
	print "$chr\n";
	
	my (%posh,%out_hash,%qual_h);
	
	#Go through each chromosome file unless it is a single chromosome dataset
	if ($species ne "single")
	{
		#$chrTemp = $chrTemp."_".$chr.".sam";
		open (IN2, "$path2output"."_$chr".".sam" ) or die "Can't open $path2output"."_$chr".".sam for reading";
	}
	else
	{
		open (IN2, "$path2sam" ) or die "Can't open $path2sam for reading";
	}
	while (my $line = <IN2>)
	{
		my @elems = split/\t/,$line;
		my $read_id = $elems[0];
		$read_id =~s/\/3$//;
		my $chrom = $elems[2];
		my $pos = $elems[3];
		my $flags = $elems[1];
		my $qual = $elems[4];
		my $strand;

		if ($qual < 10)
		{
			$qual = "0".$qual;
		}

		if (exists $posh{$read_id})
		{
			if ($qual >= 10 || $qual_h{$read_id} >= 10)
			{
				my $end = $pos + $read_length - 1;
				if ($end - $posh{$read_id} > 600)
				{
					print "Error $read_id\n";
					next;
				}
				if ($flags == 99 || $flags == 147)
				{
					$strand = "+";
				}
				elsif ($flags == 83 || $flags == 163)
				{
					$strand = "-";
				}
				else
				{
					print "Flag Problem $read_id\n";
					next;
				}
				my $string;
				if ($bed == 1)
				{
					$string = "$chrom\t$posh{$read_id}\t$end\t$read_id\t$qual_h{$read_id}$qual\t$strand\n";
				}
				else
				{
					$string = "$chrom\tSAM2pippy\tpaired\t$posh{$read_id}\t$end\t$qual_h{$read_id}$qual\t$strand\t0\t$read_id\n";
				}
				my $coords = "$chrom"."_$posh{$read_id}"."_$end";
				$count_reads++;
				my $length = $end - $posh{$read_id};
				if (exists $out_hash{$coords})
				{
					$count_rep++;
				}
				else
				{
					$out_hash{$coords} = $string;
				}
			}
			else
			{
				$failed_qc++;
			}
		}
		else
		{
			if ($se == 1)
			{
				my $end = $pos + $read_length - 1;
				if ($flags == 0)
				{
					$strand = "+";
				}
				elsif ($flags == 16)
				{
					$strand = "-";
				}
				else
				{
					$strand = ".";
				}
				my $string;
				if ($bed == 1)
				{
					$string = "$chrom\t$pos\t$end\t$read_id\t$qual\t$strand\n";
				}
				else
				{
					$string = "$chrom\tSAM2pippy\tse\t$pos\t$end\t$qual\t$strand\t0\t$read_id\n";
				}
				my $coords = "$chrom"."_$pos"."_$end";
				$count_reads++;
				if (exists $out_hash{$coords})
				{
					$count_rep++;
				}
				else
				{
					$out_hash{$coords} = $string;
				}
	
			}
			else
			{
				$posh{$read_id} = $pos;
				$qual_h{$read_id} = $qual;
			}
		}
		$count_lines ++;
	}
	close IN2;
	#$chrTemp =~s/\_$chr.sam//;
	my %new_hash;
	foreach my $key (keys %out_hash)
	{
		my @tmp = split /_/,$key;
		$new_hash{$out_hash{$key}} = $tmp[1];
	}
	foreach my $key (sort { $new_hash{$a}<=>$new_hash{$b}} keys %new_hash )
	{
		print OUT "$key";
	}
	my $percent_repeat = "NA";
	if ($count_reads > 0)
	{
		$percent_repeat = ($count_rep/$count_reads)*100;
		$percent_repeat = sprintf("%.0f", $percent_repeat);
	}
	print "$chr $count_rep $count_reads $percent_repeat%\n";

	$total_rep +=$count_rep;
	$total_read +=$count_reads;
}
my $aligned = $count_lines/2;
my $final_reads = $total_read-$total_rep;
close OUT;

my $mean_fragment = frag_length($path2output);
open (SUM ,">$path2output"."_summaryCounts.txt") or die "Can't open $path2output"."_summaryCounts.txt for writing";
print SUM "\n";
print SUM "Aligned\t$aligned\n";
print SUM "Failed QC (q<10)\t$failed_qc\n";
print SUM "High Quality aligned reads\t$total_read\n";
print SUM "Total repeated\t$total_rep\n";
print SUM "Final Reads (nonclonal)\t$final_reads\n";
print SUM "Mean Fragment length\t$mean_fragment\n";
close SUM;


foreach my $chr (@chroms)
{
	#$chrTemp = $chrTemp."_".$chr.".sam";
	my $chrTemp = "$path2output"."_$chr".".sam";
	print "unlinking $chrTemp\n";
	unlink $chrTemp or die "Can't delete $chrTemp";
	#$chrTemp =~ s/\_$chr.sam//;
}

sub frag_length
{
	my $path2file = shift;	my @fragSizeDistr;
	$fragSizeDistr[0]={}; # initialize frag size distribution hash for ith file
	my $hashRef=$fragSizeDistr[0];
	open FILE, "$path2file"; 
	my $n=0;
	my $mean_fragment = 0;
	while(<FILE>)
	{
		chop; $n++;
		my ($chr, $start, $end, @rest) = split /\t/;
		$$hashRef{$end-$start} += 1; # update frag size distribution
	}
	close FILE;
	foreach my $insSize (keys %$hashRef) {
		$$hashRef{$insSize} /= $n;
	}
	my %fragSize = %$hashRef;
	foreach my $insSize (sort numerically keys %fragSize)
	{
		$mean_fragment += $insSize * ${fragSize{$insSize}};
	}
	return $mean_fragment;
}	
sub numerically {$a<=>$b};
