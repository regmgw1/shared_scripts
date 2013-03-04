#!/usr/bin/perl -w

=head2 Authors

=head3 Written by

Reiner Schulz & Gareth Wilson
gareth.wilson@cancer.ucl.ac.uk

=head2 Description
part 2 of medusa pipeline
fragment length normalisation script obtained from Reiner Schulz, edited to reduce memory load by Gareth Wilson.
this part of the pipeline can be skipped, as long as the config file is configured correctly. The script equalises distributions between samples, however in doing so
it removes fragments from further analysis. Therefore whilst potentially improving the reliability of your dmr calls, it comes at a cost of analysing less data.

Script not fully debugged yet though distrbutions generated are normalised correctly.
=head2 Usage

Usage: ./medusa_pipe_pt2.pl norm_prefix path2beddir bedfiles (list bed file names with space inbetween)

=cut

use POSIX;
use List::Util 'shuffle';
use IO::File;
use strict;

unless (@ARGV >=3) {
        die "\n\nUsage:\n ./medusa_pipe_pt2.pl norm_prefix path2beddir bedfiles (list bed file names with space inbetween)\nPlease try again.\n\n\n";}
        
my $prefix = $ARGV[0];
my $path2files = $ARGV[1];
my @files = @ARGV[2..$#ARGV];
my (@fragSizeDistr,%minFragSizeDistr);
for( my $i=0; $i<@files; $i++) 
{
	$fragSizeDistr[$i]={}; # initialize frag size distribution hash for ith file
	my $hashRef=$fragSizeDistr[$i];
	open FILE, "$path2files/$files[$i]"; 
	my $n=0;
	while(<FILE>) 
	{
		chop; $n++;
		my ($chr, $start, $end, @rest) = split /\t/;
		$$hashRef{$end-$start} += 1; # update frag size distribution
	}
	close FILE;
	# work with relative distributions (wastes fewer reads if absolute read numbers are significantly different)
	foreach my $insSize (keys %$hashRef)
	{
		$$hashRef{$insSize} /= $n;
	}
	print STDERR "read $n records in file ${files[$i]}\n";
	if( $i == 0) 
	{ # initialize min frag size distribution
		%minFragSizeDistr=%$hashRef;
	} 
	else 
	{ # update min frag size distribution
		foreach my $insSize (keys %minFragSizeDistr) 
		{
			if( not exists $$hashRef{$insSize}) 
			{
				$minFragSizeDistr{$insSize} = 0;
			} 
			elsif( $$hashRef{$insSize} < $minFragSizeDistr{$insSize}) 
			{
				$minFragSizeDistr{$insSize} = $$hashRef{$insSize};
			}
		}
	}
}

print STDERR "min frag size distr\n";
foreach my $insSize (sort numerically keys %minFragSizeDistr) {
	print STDERR "$insSize\t${minFragSizeDistr{$insSize}}\n";
}
print STDERR "cumulative min frag size distr\n";
my $c = 0;
foreach my $insSize (sort numerically keys %minFragSizeDistr)
{
	$c+=$minFragSizeDistr{$insSize};
	print STDERR "$insSize\t$c\n";
}
# output random subsets of input BED files, each w/ the same fraction of frags for a given frag size, ie, all having same relative frag size distribution
open COUNTS, ">$path2files/fragment_counts.txt";
for(my $i=0; $i<@files; $i++) {

	my $time = time();
	print STDERR "time = $time\n";
	my %fht;
	my @tmpFileUnlink;
	foreach my $insSize (sort numerically keys %minFragSizeDistr)
	{
		my $fh = new IO::File(">$path2files/$time"."_$insSize") or die "Can't open $path2files/$time"."_$insSize";
		push @tmpFileUnlink, "$path2files/$time"."_$insSize";
		$fht{$insSize} = $fh;
	}
	#reads through file and places each fragment into a bin dependent on fragment length
	my %size2BED=();
	open FILE, "$path2files/$files[$i]"; 
	my $n=0;
	while(<FILE>) {
		chop; 
		$n++;
		my ($chr, $start, $end, @rest) = split /\t/;
		die unless $chr =~ /^chr/;
		my $insSize=$end-$start;
		if (exists $fht{$insSize})
		{
			my $fh = $fht{$insSize};
			print $fh "$chr\t$start\t$end\t@rest\n";
		}
		else
		{
			print STDERR "BUGGY: $insSize\n";
		}
	}
	close FILE;
	print STDERR "read $n records in file ${files[$i]}\n";
	foreach my $fh (values (%fht))
	{
		close $fh;
	}
	# foreach insert size, determines number of reads required, shuffles reads in that size bin and selects the relevant number of reads
	my %chr2BED=();
	open OFILE, ">$path2files/${prefix}_${files[$i]}";
	my $j = 0;
	foreach my $insSize (keys %minFragSizeDistr) {
		
		my $k= floor( $n * $minFragSizeDistr{$insSize}) -1;
		next unless $k > 0;
		open TMP, "$path2files/$time"."_$insSize" or die "Can't open $path2files/$time"."_$insSize";
		my @BEDs=<TMP>;
		close TMP;
		#unlink("$path2files/$time"."_$insSize");
		my @randBEDs=shuffle @BEDs;
		foreach my $BED (@randBEDs[0..$k]) {
			my ($chr, $start, $end, @rest) = split/\t/, $BED;
			print OFILE "$chr\t$start\t$end\t@rest";
			$j++;
		}
	}
	unlink(@tmpFileUnlink);
	# opens output file, then sorts the data according to chr, start, stop etc.
	
	close OFILE;
	
	print COUNTS "${prefix}_${files[$i]}\t$n\n";
	print STDERR "wrote $j records in file ${prefix}_${files[$i]}\n";
	
}
close COUNTS;

sub numerically {$a<=>$b};
