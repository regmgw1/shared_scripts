#!/usr/bin/perl

# fragment length normalisation script obtained from Reiner Schulz, minor edits by me.

use POSIX;
use List::Util 'shuffle';

$prefix = $ARGV[0];
$path2files = $ARGV[1];
@files = @ARGV[2..$#ARGV];
#print STDERR join( ",", @files), "\n";
for( $i=0; $i<@files; $i++) {
	$fragSizeDistr[$i]={}; # initialize frag size distribution hash for ith file
	my $hashRef=$fragSizeDistr[$i];
	open FILE, "$path2files/$files[$i]"; $n=0;
	while(<FILE>) {
		chop; $n++;
		($chr, $start, $end, @rest) = split /\t/;
		$$hashRef{$end-$start} += 1; # update frag size distribution
	}
	close FILE;
	# work with relative distributions (wastes fewer reads if absolute read numbers are significantly different)
	foreach $insSize (keys %$hashRef) {
		$$hashRef{$insSize} /= $n;
	}
	print STDERR "read $n records in file ${files[$i]}\n";
	if( $i == 0) { # initialize min frag size distribution
		%minFragSizeDistr=%$hashRef;
	} else { # update min frag size distribution
		foreach $insSize (keys %minFragSizeDistr) {
			if( not exists $$hashRef{$insSize}) {
				$minFragSizeDistr{$insSize} = 0;
			} elsif( $$hashRef{$insSize} < $minFragSizeDistr{$insSize}) {
				$minFragSizeDistr{$insSize} = $$hashRef{$insSize};
			}
		}
	}
}
print STDERR "min frag size distr\n";
foreach $insSize (sort numerically keys %minFragSizeDistr) {
	print STDERR "$insSize\t${minFragSizeDistr{$insSize}}\n";
}
print STDERR "cumulative min frag size distr\n";
$c = 0;
foreach $insSize (sort numerically keys %minFragSizeDistr) {
	$c+=$minFragSizeDistr{$insSize};
	print STDERR "$insSize\t$c\n";
}
# output random subsets of input BED files, each w/ the same fraction of frags for a given frag size, ie, all having same relative frag size distribution
open COUNTS, ">$path2files/fragment_counts.txt";
for( $i=0; $i<@files; $i++) {

	my %size2BED=();
	open FILE, "$path2files/$files[$i]"; $n=0;
	while(<FILE>) {
		chop; $n++;
		($chr, $start, $end, @rest) = split /\t/;
		die unless $chr =~ /^chr/;
		$insSize=$end-$start;
		$size2BED{$insSize} = [] unless exists $size2BED{$insSize};
		my $BEDs=$size2BED{$insSize};
		push @$BEDs, [$chr, $start, $end, @rest];
	}
	close FILE;
	print STDERR "read $n records in file ${files[$i]}\n";
	
	my %chr2BED=();
	foreach $insSize (keys %minFragSizeDistr) {
		$k= floor( $n * $minFragSizeDistr{$insSize}) -1;
		next unless $k > 0;
		my $BEDs=$size2BED{$insSize};
		my @randBEDs=shuffle @$BEDs;
#		print STDERR scalar @randBEDs, "\n";
		foreach my $BED (@randBEDs[0..$k]) {
			($chr, $start, $end, @rest) = @$BED;
			$chr2BED{$chr}= {} unless exists $chr2BED{$chr};
			my $start2BED= $chr2BED{$chr};
			$$start2BED{$start}= {} unless exists $$start2BED{$start};
			my $end2BED= $$start2BED{$start};
			$$end2BED{$end}= [] unless exists $$end2BED{$end};
			my $BEDs= $$end2BED{$end};
			push @$BEDs, $BED;
		}
	}
	open OFILE, ">$path2files/${prefix}_${files[$i]}"; $n=0;
	foreach $chr (sort keys %chr2BED) {
		my $start2BED= $chr2BED{$chr};
		foreach $start (sort numerically keys %$start2BED) {
			my $end2BED= $$start2BED{$start};
			foreach $end (sort numerically keys %$end2BED) {
				my $BEDs= $$end2BED{$end};
				foreach $BED (@$BEDs) {
					($chr, $start, $end, @rest) = @$BED;
					print OFILE join( "\t", @$BED), "\n"; $n++;
				}
			}
		}
	}
	close OFILE;
	
	print COUNTS "${prefix}_${files[$i]}\t$n\n";
	print STDERR "wrote $n records in file ${prefix}_${files[$i]}\n";
}
close COUNTS;
sub numerically {$a<=>$b};
