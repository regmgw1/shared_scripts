#!/bin/bash

USEQ_HOME=/usr/local/bioinf/USeq_6.8/Apps

# $1 = comma-separated list of treatment point data folders
# $2 = comma-separated list of control point data folders
# $3 = output folder name (must not exist)
# $4 = input directory
# $5 = minimum_read_depth (if in doubt use default of 10)
# $6 = minimum FDR threshold for saving windows (if in doubt use default 0.5)

CURRENT=`pwd`
cd $4

java -Xmx10G -jar $USEQ_HOME/MultipleReplicaScanSeqs -t $1 -c $2 -s $3 -p 0 -w 500 -r /usr/local/bin/R -m $5 -f $6
# call DMRs (not using sub-window option of EnrichedRegionMaker: not relevant for MeDIP-seq)
# hypermethylation
# 20% FDR
java -Xmx2G -jar $USEQ_HOME/EnrichedRegionMaker -f $3 -s 7,0.001 -i 0,1
# 10% FDR
java -Xmx2G -jar $USEQ_HOME/EnrichedRegionMaker -f $3 -s 10,0.001 -i 0,1
# 5% FDR
java -Xmx2G -jar $USEQ_HOME/EnrichedRegionMaker -f $3 -s 13,0.001 -i 0,1
# 1% FDR
java -Xmx2G -jar $USEQ_HOME/EnrichedRegionMaker -f $3 -s 20,0.001 -i 0,1
# hypomethylation
# 20% FDR
java -Xmx2G -jar $USEQ_HOME/EnrichedRegionMaker -f $3 -s 7,0.001 -i 0,1 -m
# 10% FDR
java -Xmx2G -jar $USEQ_HOME/EnrichedRegionMaker -f $3 -s 10,0.001 -i 0,1 -m
# 5% FDR
java -Xmx2G -jar $USEQ_HOME/EnrichedRegionMaker -f $3 -s 13,0.001 -i 0,1 -m
# 1% FDR
java -Xmx2G -jar $USEQ_HOME/EnrichedRegionMaker -f $3 -s 20,0.001 -i 0,1 -m


# TODO: extract report ... something like this:
echo "hypERmethylated DMRs:"
for fdr in 7 10 13 20; do ls -d $3/Enr*FDR$fdr.0_Log*/ |sed "s/\([^_]*\)_MRSS\/EnrichedRegions_BinaryWindowData_FDR$fdr\.0_Log2Ratio0\.001_\([0-9]*\)\//\1\t\2/"; done |perl -e 'while(<>){chop;($s,$n)=split /\s+/;$s2n{$s}=[] unless exists $s2n{$s};$aref=$s2n{$s};push @$aref, $n;};foreach $s (sort keys %s2n){$aref=$s2n{$s};print "$s\t",join("\t",@$aref),"\n";};'
echo "hypOmethylated DMRs:"
for fdr in 7 10 13 20; do ls -d $3/Red*FDR$fdr.0_Log*/ |sed "s/\([^_]*\)_MRSS\/ReducedRegions_BinaryWindowData_FDR$fdr\.0_Log2Ratio0\.001_\([0-9]*\)\//\1\t\2/"; done |perl -e 'while(<>){chop;($s,$n)=split /\s+/;$s2n{$s}=[] unless exists $s2n{$s};$aref=$s2n{$s};push @$aref, $n;};foreach $s (sort keys %s2n){$aref=$s2n{$s};print "$s\t",join("\t",@$aref),"\n";};'

cd $CURRENT
