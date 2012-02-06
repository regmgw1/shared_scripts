#!/bin/bash

# part 0 of medusa pipeline - written by Gareth Wilson
# runs initial alignment and filtering. Currently performed using BWA (http://bio-bwa.sourceforge.net/) and samtools (http://samtools.sourceforge.net/).
# limited alignment parameters are passed from medusa.pl. defaults have generally worked well for our data. More parameters can be added if required (edit PARAM_BWA_ALN or PARAM_BWA_SAMPE).
# in future no reason why other alignment programme can't be used - just need to ensure output is in sam format.

SAMPLE=$1
PATH2READS=$2
PATH2GENOME=$3
READ1=$4
READ2=$5
BWAGENOME=$6
SAMGENOME=$7
PATH2OUTPUT=$8
MAX_INSERT=$9

INPUT1=$PATH2READS/$READ1
INPUT2=$PATH2READS/$READ2
IN_GENOME=$PATH2GENOME/$BWAGENOME
IN_S_GENOME=$PATH2GENOME/$SAMGENOME

PARAM_BWA_ALN="$IN_GENOME" 
PARAM_BWA_SAMPE="-a $MAX_INSERT -P $IN_GENOME"
PARAM_SAM="-S -t $IN_S_GENOME"".fai -b"

bwa aln $PARAM_BWA_ALN $INPUT1 >$PATH2OUTPUT/$SAMPLE""_1.sai
bwa aln $PARAM_BWA_ALN $INPUT2 >$PATH2OUTPUT/$SAMPLE""_2.sai
bwa sampe $PARAM_BWA_SAMPE $PATH2OUTPUT/$SAMPLE""_1.sai $PATH2OUTPUT/$SAMPLE""_2.sai $INPUT1 $INPUT2 >$PATH2OUTPUT/$SAMPLE"".sam
samtools view $PARAM_SAM $PATH2OUTPUT/$SAMPLE"".sam >$PATH2OUTPUT/$SAMPLE"".bam
samtools sort $PATH2OUTPUT/$SAMPLE"".bam $PATH2OUTPUT/$SAMPLE"".sorted
samtools index $PATH2OUTPUT/$SAMPLE"".sorted.bam
samtools view -f 2 $PATH2OUTPUT/$SAMPLE"".sorted.bam >$PATH2OUTPUT/$SAMPLE""_filter.sam

