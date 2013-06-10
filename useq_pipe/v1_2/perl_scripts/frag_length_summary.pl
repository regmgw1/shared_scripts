#!/usr/bin/perl -w

=head2 Authors

=head3 Written by

Reiner Schulz & Gareth Wilson
gareth.wilson@cancer.ucl.ac.uk

=head2 Description
fragment length normalisation script originally written by Reiner Schulz, gross edits by Gareth Wilson have reduced it to a script to generate fragment length distributions.
not currently part of medusa, but would be useful to add.
=head2 Usage

Usage: ./frag_length_summary.pl path2genome path2peakroot path2peaklist path2output

=cut

# 

use POSIX;
use List::Util 'shuffle';
use strict;

unless (@ARGV ==1) {
        die "\n\nUsage:\n ./frag_length_summary.pl path2bedFile\nPlease try again.\n\n\n";}

my $path2file = $ARGV[0];open (OUT, ">$path2file"."_frag_dist.txt" ) or die "Can't open $path2file"."_frag_dist.txt for writing";
my @fragSizeDistr;
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
print STDERR "read $n records in file $path2file\n";
my %fragSize = %$hashRef;
foreach my $insSize (sort numerically keys %fragSize)
{
	print OUT "$insSize\t${fragSize{$insSize}}\n";
	$mean_fragment += $insSize * ${fragSize{$insSize}};
}
close OUT;
print "Mean Fragment Length $path2file = $mean_fragment\n";
sub numerically {$a<=>$b};
