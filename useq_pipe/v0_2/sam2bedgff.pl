#!/usr/bin/perl -w


#################################################################
# sam2bedgff.pl - to be run on filtered sam files. convert to bed/gff with only unique reads included
#################################################################

use strict;
use File::Basename;

unless (@ARGV ==7 ) {
        die "\n\nUsage:\n ./sam2pippy_v3.pl path2sam sample species read_length single_end? bed? path2output\nPlease try again.\n\n\n";}

my $path2sam = shift;
my $sample = shift;
my $species = shift;
my $read_length = shift;
my $se = shift;
my $bed = shift;
my $path2output = shift;

my @chroms;
if ($species eq "Mouse")
{
	@chroms = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,'X','Y');
}
elsif ($species eq "Human")
{
	@chroms = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y');
}
elsif ($species eq "single")
{
	@chroms = (1);
}
else
{
	die "Currently only configured for Mouse and Human";
}
my $total_rep = 0;
my $total_read = 0;
my $failed_qc = 0;

if ($bed == 1)
{
	open (OUT ,">$path2output") or die "Can't open $path2output for writing";
}
else
{
	open (OUT ,">$path2output") or die "Can't open $path2output for writing";
}

foreach my $chr (@chroms)
{
	my $count_rep = 0;
	my $count_reads = 0;
	print "$chr\n";
	
	my (%posh,%out_hash,%qual_h);
	open (IN, "$path2sam" ) or die "Can't open $path2sam for reading";
	while (my $line = <IN>)
	{
		my @elems = split/\t/,$line;
		my $read_id = $elems[0];
		my $chrom = $elems[2];
		my $pos = $elems[3];
		my $flags = $elems[1];
		my $qual = $elems[4];
		my $strand;
		if ($qual < 10)
		{
			$qual = "0".$qual;
		}
		if ($chrom eq $chr || $species eq "single")
		{
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
					$strand = ".";
				}
				my $string;
				if ($bed == 1)
				{
					$string = "chr$chrom\t$posh{$read_id}\t$end\t$read_id\t$qual_h{$read_id}$qual\t$strand\n";
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
		}
		else
		{
			next;
		}
	}
	close IN;
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
	print "$chr $count_rep $count_reads\n";
	
	$total_rep +=$count_rep;
	$total_read +=$count_reads;
}
print "Total reads = $total_read\n";
print "Total repeated = $total_rep\n";
print "Failed QC (q<10) = $failed_qc\n";
close OUT;
