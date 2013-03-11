#!/usr/bin/perl -w

=head2 Authors

=head3 Created by

Gareth Wilson
gareth.wilson@cancer.ucl.ac.uk

=head2 Description
Use bedTools 'fastaFromBed' to obtain metadata for the useq peaks
=head2 Usage

Usage: ./peaks_cg_bedtools.pl path2genome path2peakroot path2peaklist path2output

=cut

#################################################################
# peaks_cg_bedtools.pl
#################################################################

use strict;

unless (@ARGV ==7) {
        die "\n\nUsage:\n ./peaks_cg_bedtools.pl path2genome path2peakroot path2peaklist path2genes upstream_threshold downstream_threshold path2output\nPlease try again.\n\n\n";}

my $path2genome = shift;
my $path2peakroot = shift;
my $path2peaklist = shift;
my $path2genes = shift;
my $up = shift;
my $down = shift;
my $path2output = shift;

my (@samples, @peak_subs);

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
	open (OUT, ">$path2output/$peaksub"."_peak_cg.txt") or die "Can't open $path2output/$peaksub"."_peak_cg.txt for writing";
	print OUT "DMR_chr\tDMR_start\tDMR_stop\tLength\tCpG_count\tCpG_density\tC_count\tC_density\tG_count\tG_density\tBH_Cords\tBH_Length\tBH_CpG_count\tBH_CpG_density\tBH_C_count\tBH_C_density\tBH_G_count\tBH_G_density\tBH_FDR\tBH_Log2\tNear_gene_ID_up_$up"."_down_$down\tGene_coords\tStrand\n";
	my $peakfile = $peaksub;
	my (%hash,%bh_hash, %scores);
	$peakfile =~s/.*Binary/binary/;
	open (TMP, ">$path2output/tmp.gff") or die "Can't open $path2output/tmp.gff for writing";
	open (BHTMP, ">$path2output/tmpBH.gff") or die "Can't open $path2output/tmpBH.gff for writing";
	my $l_count = 0;
	my $nearest = nearest_gene("$path2peakroot/$peaksub/$peakfile".".gff",$path2genes,$up,$down); 
	my %nearHash = %$nearest;
	open (IN, "$path2peakroot/$peaksub/$peakfile".".xls" ) or die "Can't open $path2peakroot/$peaksub/$peakfile".".xls for reading";
	while (my $line = <IN>)
	{
		if ($l_count > 0)
		{
			chomp $line;
			my @elems = split/\t/, $line;
			my $tmpC = $elems[1];
			$tmpC =~s/^chr//;
			print TMP "chr$tmpC\ttmpGFF\tmoreTemp\t$elems[2]\t$elems[3]\t$elems[4]\t.\t.\tnoMore\n";
			print BHTMP "chr$tmpC\ttmpGFF\tmoreTemp\t$elems[5]\t$elems[6]\t$elems[7]\t$elems[8]\t.\tnoMore\n";
			my $start = $elems[2] - 1;
			my $startBH = $elems[5] - 1;
			my $coords = "chr$tmpC".":$start"."-$elems[3]";
			my $coordsBH = "chr$tmpC".":$startBH"."-$elems[6]";
			$hash{$coordsBH} = $coords;
			$scores{$coordsBH} = "$elems[7]\t$elems[8]";
		}
		$l_count++;
	}
	close IN;
	close TMP;
	close BHTMP;
	my @seq = `fastaFromBed -fi $path2genome -bed  $path2output/tmpBH.gff -fo $path2output/tmpBH.fasta -tab`;
	open (IN, "$path2output/tmpBH.fasta" ) or die "Can't open $path2output/tmpBH.fasta for reading";
	while (my $line = <IN>)
	{
		chomp $line;
		my @elems = split/\t/,$line;
		my $seq = $elems[1];
		my @tmp = ($seq =~/CG/g);
		my $cgs = $#tmp + 1;
		my $seq_length = length($seq);
		my $density = ($cgs/$seq_length * 2) * 100;
		@tmp = ($seq =~/C/g);
		my $cs = $#tmp + 1;
		my $c_density = ($cs/$seq_length) * 100;
		@tmp = ($seq =~/G/g);
		my $gs = $#tmp + 1;
		my $g_density = ($gs/$seq_length) * 100;
		my $id_whole = $hash{$elems[0]};
		my $scores = $scores{$elems[0]};
		$bh_hash{$id_whole} = "$elems[0]\t$seq_length\t$cgs\t$density\t$cs\t$c_density\t$gs\t$g_density\t$scores";
	}
	close IN;
	my @seq2 = `fastaFromBed -fi $path2genome -bed  $path2output/tmp.gff -fo $path2output/tmp.fasta -tab`;
	open (IN, "$path2output/tmp.fasta" ) or die "Can't open $path2output/tmp.fasta for reading";
	while (my $line = <IN>)
	{
		chomp $line;
		my @elems = split/\t/,$line;
		my $seq = $elems[1];
		my @tmp = ($seq =~/CG/g);
		my $cgs = $#tmp + 1;
		my $seq_length = length($seq);
		my $density = ($cgs/$seq_length * 2) * 100;
		@tmp = ($seq =~/C/g);
		my $cs = $#tmp + 1;
		my $c_density = ($cs/$seq_length) * 100;
		@tmp = ($seq =~/G/g);
		my $gs = $#tmp + 1;
		my $g_density = ($gs/$seq_length) * 100;
		my $bt_coords = $elems[0];
		my ($outChr, $outStart, $outStop);
		if ($bt_coords =~m/(\w+):(\d+)-(\d+)/)
		{
			$outChr = $1;
			$outStart = $2 + 1;
			$outStop = $3;
		}
		else
		{
			die "problem with coord regex\n";
		}
		if (exists $nearHash{$elems[0]})
		{
			print OUT "$outChr\t$outStart\t$outStop\t$seq_length\t$cgs\t$density\t$cs\t$c_density\t$gs\t$g_density\t$bh_hash{$elems[0]}\t$nearHash{$elems[0]}\n";
		}
		else
		{
			print "elems = $elems[0]\n";			
			print OUT "$outChr\t$outStart\t$outStop\t$seq_length\t$cgs\t$density\t$cs\t$c_density\t$gs\t$g_density\t$bh_hash{$elems[0]}\t\t\t\n";
		}
	}
	close IN;
	close OUT;
	unlink ("$path2output/tmp.fasta","$path2output/tmp.gff","$path2output/tmpBH.fasta","$path2output/tmpBH.gff","$path2output/tmpNear.gff");
}

sub nearest_gene
{
	my $path2peaks = shift;
	my $path2genes = shift;
	my $up = shift;
	my $down = shift;
	
	my %hash;
	my %inter;
	my %inter_dup;
	my %nearHash;

	open (GFF, ">$path2output/tmpNear.gff") or die "Can't open $path2output/tmpNear.gff for writing";
	open (DMR, $path2peaks ) or die "Can't open $path2peaks for reading";
	while (my $dmr = <DMR>)
	{
		$dmr =~s/^chr//;		
		print GFF "$dmr";		
	}
	close GFF;
	my @intersect = `intersectBed -a $path2genes -b $path2output/tmpNear.gff -wa -wb`;
	my @windows = `windowBed -a $path2genes -b $path2output/tmpNear.gff -l $up -r $down -sw`;

	foreach my $inter (@intersect)
	{
		chomp $inter;
		my @elems = split /\t/,$inter;
		my $coords = "$elems[0]".":$elems[12]"."-$elems[13]";
		if (exists $inter{$coords})
		{
			$inter_dup{$coords} = $inter{$coords};
		}
		$inter{$coords} = $inter;
	}
	

	foreach my $win (@windows)
	{
		chomp $win;
		my @elems = split /\t/,$win;
		my $coords = "$elems[0]".":$elems[12]"."-$elems[13]";
		my $start = $elems[12];
		my ($new_dist, $old_dist);
		if (exists $hash{$coords})
		{
			if ($elems[6] eq "+")
			{
				$new_dist = abs($elems[3] - $start);
			}
			else
			{
				$new_dist = abs($elems[4] - $start);
			}
			my @old_elems = split/\t/,$hash{$coords};
			if ($old_elems[6] eq "+")
			{
				$old_dist = abs($old_elems[3] - $start);
			}
			else
			{
				$old_dist = abs($old_elems[4] - $start);
			}
			if ($new_dist < $old_dist)
			{
				$hash{$coords} = $win;
			}
		}
		else
		{
			$hash{$coords} = $win;
		}
		if (exists $inter{$coords})
		{
			$hash{$coords} = $inter{$coords};
		}
	}
	foreach my $outdup (keys %inter_dup)
	{
		my @elems = split/\t/,$inter_dup{$outdup};
		my $mod_start = $elems[12] - 1;
		my $coords = "$elems[0]".":$mod_start"."-$elems[13]";
		my @waste = split/_/,$elems[1];
		my $id = $waste[$#waste];
		$nearHash{$coords} = "$id\t$elems[2]\t$elems[6]";
	}
	foreach my $out (keys %hash)
	{
		my @elems = split/\t/,$hash{$out};
		my $mod_start = $elems[12] - 1;
		my $coords = "chr$elems[0]".":$mod_start"."-$elems[13]";
		my @waste = split/_/,$elems[1];
		my $id = $waste[$#waste];
		$nearHash{$coords} = "$id\t$elems[2]\t$elems[6]";
		print "$coords\t$nearHash{$coords}\n";
	}
	return \%nearHash;
}	
