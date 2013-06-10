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

unless (@ARGV ==8) {
        die "\n\nUsage:\n ./peaks_cg_bedtools.pl path2genome path2peakroot path2peaklist path2genes path2chromlist upstream_threshold downstream_threshold path2output\nPlease try again.\n\n\n";}

my $path2genome = shift;
my $path2dmrs = shift;
my $path2peaklist = shift;
my $path2genes = shift;
my $path2chrom = shift;
my $up = shift;
my $down = shift;
my $path2output = shift;

my (@samples, @dmr_files);

#open chrom list, determine nomenclature for sequence retrieval
my $chrom_nom = "";
open (CHR, "$path2chrom" ) or die "Can't open $path2chrom for reading";
while (my $line = <CHR>)
{
	chomp $line;	
	if ($line =~m/^chr/)
	{
		$chrom_nom = "chr";
	}	
}
close CHR;

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
	open (OUT, ">$path2output/$dmr_file"."_peak_cg.txt") or die "Can't open $path2output/$dmr_file"."_peak_cg.txt for writing";
	print OUT "DMR_chr\tDMR_start\tDMR_stop\tLength\tCpG_count\tCpG_density\tC_count\tC_density\tG_count\tG_density\tAdj.pvalue\tLogFC\tNear_gene_ID_up_$up"."_down_$down\tGene_coords\tStrand\n";
	my (%hash, %scores);
	open (TMP, ">$path2output/tmp.gff") or die "Can't open $path2output/tmp.gff for writing";
	my $l_count = 0;
	my $nearest = nearest_gene("$path2dmrs/$dmr_file.txt",$path2genes,$up,$down,$chrom_nom); 
	my %nearHash = %$nearest;
	open (IN, "$path2dmrs/$dmr_file.txt" ) or die "Can't open $path2dmrs/$dmr_file.txt for reading";
	while (my $line = <IN>)
	{
		if ($l_count > 0)
		{
			chomp $line;
			my @elems = split/\t/, $line;
			my $tmpC = $elems[0];
			$tmpC =~s/^chr//;
			print TMP "$chrom_nom$tmpC\ttmpGFF\tmoreTemp\t$elems[1]\t$elems[2]\t$elems[3]\t.\t.\tnoMore\n";
			my $start = $elems[1] - 1;
			my $coords = "$chrom_nom$tmpC".":$start"."-$elems[2]";
			$hash{$coords} = $coords;
			$scores{$coords} = "$elems[31]\t$elems[28]";
		}
		$l_count++;
	}
	close IN;
	close TMP;
	my @seq = `fastaFromBed -fi $path2genome -bed  $path2output/tmp.gff -fo $path2output/tmp.fasta -tab`;
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
		if ($bt_coords =~m/(\S+):(\d+)-(\d+)/)
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
			print OUT "$outChr\t$outStart\t$outStop\t$seq_length\t$cgs\t$density\t$cs\t$c_density\t$gs\t$g_density\t$scores{$elems[0]}\t$nearHash{$elems[0]}\n";
		}
		else
		{
			print OUT "$outChr\t$outStart\t$outStop\t$seq_length\t$cgs\t$density\t$cs\t$c_density\t$gs\t$g_density\t$scores{$elems[0]}\t\t\t\n";
		}
	}
	close IN;
	close OUT;
	unlink ("$path2output/tmp.fasta","$path2output/tmp.gff","$path2output/tmpBH.fasta","$path2output/tmpBH.gff","$path2output/tmpNear.bed");
}
sub nearest_gene
{
	my $path2peaks = shift;
	my $path2genes = shift;
	my $up = shift;
	my $down = shift;
	my $chromNom = shift;
	my %hash;
	my %inter;
	my %inter_dup;
	my %nearHash;
	my $dmrC = 0;
	open (BED, ">$path2output/tmpNear.bed") or die "Can't open $path2output/tmpNear.bed for writing";
	open (DMR, $path2peaks ) or die "Can't open $path2peaks for reading";
	while (my $dmr = <DMR>)
	{
		if ($dmrC > 0)
		{
			$dmr =~s/^chr//;
			my @elems = split/\t/,$dmr;		
			print BED "$elems[0]\t$elems[1]\t$elems[2]\t$elems[3]\n";
		}
		$dmrC++;		
	}
	close BED;
	my @intersect = `intersectBed -a $path2genes -b $path2output/tmpNear.bed -wa -wb`;
	my @windows = `windowBed -a $path2genes -b $path2output/tmpNear.bed -l $up -r $down -sw`;
	foreach my $inter (@intersect)
	{
		chomp $inter;
		my @elems = split /\t/,$inter;
		my $coords = "$elems[0]".":$elems[10]"."-$elems[11]";
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
		my $coords = "$elems[0]".":$elems[10]"."-$elems[11]";
		my $start = $elems[10];
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
		my $mod_start = $elems[10] - 1;
		my $coords = "$chromNom$elems[0]".":$mod_start"."-$elems[11]";
		my @waste = split/_/,$elems[1];
		my $id = $waste[$#waste];
		$nearHash{$coords} = "$id\t$elems[2]\t$elems[6]";
	}
	foreach my $out (keys %hash)
	{
		my @elems = split/\t/,$hash{$out};
		my $mod_start = $elems[10] - 1;
		my $coords = "$chromNom$elems[0]".":$mod_start"."-$elems[11]";
		my @waste = split/_/,$elems[1];
		my $id = $waste[$#waste];
		$nearHash{$coords} = "$id\t$elems[2]\t$elems[6]";
	}
	return \%nearHash;
}	
