#!/usr/bin/perl -w


#################################################################
# sam2bedgff.pl - to be run on filtered sam files. convert to bed/gff with only unique reads included
#################################################################


#command used to run: perl scriptsGareth/sam2bedgff.pl newTest/sample1_filter.sam sample1 Mouse 36 0 1 newTest/Output/sample1.bed
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

#Determining number of chromosomes based on species
my @chroms;
if ($species eq "Mouse")
{
	@chroms = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,'X','Y');
}
elsif ($species eq "Human")
{
	@chroms = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y');
}
elsif ($species eq "Dog")
{
	@chroms = (1..38,'X',);
}
elsif ($species eq "Chimp")
{
	@chroms = (1,'2a','2b',3..22,'X','Y');
}
elsif ($species eq "Macaca")
{
	@chroms = (1..20,'X')
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
my $count_lines = 0;

open (OUT ,">$path2output") or die "Can't open $path2output for writing";

#Separate chromosomes into individual temporary files
my $chrTemp;
#$chrTemp = $path2sam;
#$chrTemp =~ s/\.sam/\_chr/;
$chrTemp = $path2output;
$chrTemp =~ s/\.(bed|gff)/\_chr/;

print "Creating temporary chromosome files\n";
if ($species ne "single")
{
	open (IN, "$path2sam" ) or die "Can't open $path2sam for reading";
	while (my $line = <IN>)
	{
		my @elems = split/\t/,$line;
		my $chrom = $elems[2];

		$chrTemp = $chrTemp."_".$chrom.".sam";
		open (OUT2 ,">>$chrTemp") or die "Can't open $chrTemp for writing";
		print OUT2 $line;
		close OUT2;
		$chrTemp =~ s/\_$chrom.sam//;
	}
	close IN;
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
		$chrTemp = $chrTemp."_".$chr.".sam";
		open (IN2, "$chrTemp" ) or die "Can't open $chrTemp for reading";
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
					$string = "chr$chrom\t$pos\t$end\t$read_id\t$qual\t$strand\n";
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
	$chrTemp =~ s/\_$chr.sam//;
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
	my $percent_repeat = ($count_rep/$count_reads)*100;
	$percent_repeat = sprintf("%.0f", $percent_repeat);
	print "$chr $count_rep $count_reads $percent_repeat%\n";

	$total_rep +=$count_rep;
	$total_read +=$count_reads;
}
my $aligned = $count_lines/2;
my $final_reads = $total_read-$total_rep;

print "\n"; #Total Reads = \n";
print "Aligned = $aligned\n";
print "Failed QC (q<10) = $failed_qc\n";
print "High Quality aligned reads = $total_read\n";
print "Total repeated = $total_rep\n";
print "Final Reads (nonclonal) = $final_reads\n";
close OUT;

foreach my $chr (@chroms)
{
	$chrTemp = $chrTemp."_".$chr.".sam";
	unlink $chrTemp or die "Can't delete $chrTemp";
	$chrTemp =~ s/\_$chr.sam//;
}
