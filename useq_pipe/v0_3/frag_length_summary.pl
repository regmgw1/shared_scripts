#!/usr/bin/perl

# fragment length normalisation script originally written by Reiner Schulz, gross edits by me (Gareth Wilson) have reduced it to a script to generate fragment length distributions.

use POSIX;
use List::Util 'shuffle';

$prefix = $ARGV[0];
$path2files = $ARGV[1];
@files = @ARGV[2..$#ARGV];
#print STDERR join( ",", @files), "\n";
for( $i=0; $i<@files; $i++) {
	open (OUT, ">$files[$i]"."_frag_dist.txt" ) or die "Can't open $files[$i]"."_frag_dist.txt for writing";
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
	my %fragSize = %$hashRef;
	foreach $insSize (sort numerically keys %fragSize)
	{
		print OUT "$insSize\t${fragSize{$insSize}}\n";
	}
	close OUT;
}


sub numerically {$a<=>$b};
