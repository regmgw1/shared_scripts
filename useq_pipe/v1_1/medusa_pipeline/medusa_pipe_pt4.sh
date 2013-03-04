#!/bin/bash

# part 4 of medusa pipeline - written by Gareth Wilson and Reiner Schulz
# responsible for running USeq MultipleReplicaScanSeqs and USeq EnrichedRegionMaker. Need to be aware of memory requirements. More samples and more reads = more memory!
# EnrichedRegionMaker is run at 4 FDR thresholds: 20% (7), 10% (10), 5% (13) and 1% (20).

# $1 = comma-separated list of treatment point data folders
# $2 = comma-separated list of control point data folders
# $3 = output folder name (must not exist)
# $4 = input directory
# $5 = minimum_read_depth (if in doubt use default of 10)
# $6 = minimum FDR threshold for saving windows (if in doubt use default 0.5)
# $7 = bypass_variance_outlier_filter
# $8 = path to useq installation


USEQ_HOME=$8

CURRENT=`pwd`
cd $4

# bypass the variance outlier filter
if [ "$7" == "0" ]; then
java -Xmx10G -jar $USEQ_HOME/MultipleReplicaScanSeqs -t $1 -c $2 -s $3 -p 0 -w 500 -r /usr/local/bin/R -m $5 -f $6  
elif [ "$7" == "1" ]; then
java -Xmx10G -jar $USEQ_HOME/MultipleReplicaScanSeqs -t $1 -c $2 -s $3 -p 0 -w 500 -r /usr/local/bin/R -m $5 -f $6 -b
fi

java -Xmx10G -jar $USEQ_HOME/MultipleReplicaScanSeqs -t $1 -c $2 -s $3 -p 0 -w 500 -r /usr/local/bin/R -m $5 -f $6
# call DMRs (not using sub-window option of EnrichedRegionMaker: not relevant for MeDIP-seq)
# hypermethylation
# 20% FDR
java -Xmx6G -jar $USEQ_HOME/EnrichedRegionMaker -f $3/* -s 7,0.001 -i 0,1
# 10% FDR
java -Xmx6G -jar $USEQ_HOME/EnrichedRegionMaker -f $3/* -s 10,0.001 -i 0,1
# 5% FDR
java -Xmx6G -jar $USEQ_HOME/EnrichedRegionMaker -f $3/* -s 13,0.001 -i 0,1
# 1% FDR
java -Xmx6G -jar $USEQ_HOME/EnrichedRegionMaker -f $3/* -s 20,0.001 -i 0,1
# hypomethylation
# 20% FDR
java -Xmx6G -jar $USEQ_HOME/EnrichedRegionMaker -f $3/* -s 7,0.001 -i 0,1 -m
# 10% FDR
java -Xmx6G -jar $USEQ_HOME/EnrichedRegionMaker -f $3/* -s 10,0.001 -i 0,1 -m
# 5% FDR
java -Xmx6G -jar $USEQ_HOME/EnrichedRegionMaker -f $3/* -s 13,0.001 -i 0,1 -m
# 1% FDR
java -Xmx6G -jar $USEQ_HOME/EnrichedRegionMaker -f $3/* -s 20,0.001 -i 0,1 -m

echo "hypERmethylated DMRs output directories:"
for fdr in 7 10 13 20; do ls -d $3/Enr*FDR$fdr.0_Log*/;done
echo "hypOmethylated DMRs output directories:"
for fdr in 7 10 13 20; do ls -d $3/Red*FDR$fdr.0_Log*/;done

cd $CURRENT
