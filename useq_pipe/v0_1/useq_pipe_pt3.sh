#!/bin/bash

USEQ_HOME=/usr/local/bioinf/USeq_6.8/Apps
PERL_SCRIPTS=/home/POGBDOM/regmgw1/useq_pipe/
REFNAME=mm9
FILTERED_NAME=`echo $1 |sed 's/\.sam$//'`
PATH2INPUT=$2

CURRENT=`pwd`
cd $PATH2INPUT
# convert to USeq binary mid-point format
# creates folder ${NAME}_Point
java -Xmx6G -jar $USEQ_HOME/Tag2Point -f $FILTERED_NAME.bed -v $REFNAME -b -i 5 
# turn mid-point data into UCSC bigWig track data
# sum duplicate points and merge strands
# do NOT use for MultipleReplicaScanSeqs: wants stranded point data, and we want to keep "duplicates" since they are distinct PE reads w/ the same mid-point
java -Xmx6G -jar $USEQ_HOME/PointDataManipulator -p ${FILTERED_NAME}_Point -s ${FILTERED_NAME}_MergedPoint -o -i -m
# unzip BAR files, convert to textual GR and then to variable step WIG
for ZIP in ${FILTERED_NAME}_MergedPoint/*.bar.zip
do
unzip $ZIP -d ${FILTERED_NAME}_MergedPoint/
done
java -Xmx16G -jar $USEQ_HOME/Bar2Gr -f ${FILTERED_NAME}_MergedPoint/
# concat GR files and normalize scores by number of million mapped reads (points in BAR files)
for GR in ${FILTERED_NAME}_MergedPoint/*.gr
do
echo $GR
#CHR=`echo $GR |sed 's/[^/]*\/\([^_]*\).*/chr\1/'`
CHR=`echo $GR |sed 's/[^/]*\/\(chr[^_]*\).*/\1/'`
echo $CHR
echo "variableStep chrom=$CHR" >>${FILTERED_NAME}_MergedPoint/${FILTERED_NAME}.wig
NUMPTS=`grep $FILTERED_NAME fragment_counts.txt|cut -f 2`
perl -e "while(<>){chop;(\$pos,\$s)=split /\s+/;print \"\$pos\t\",sprintf(\"%0.4f\",\$s/($NUMPTS/1000000)),\"\n\";}" <$GR >>${FILTERED_NAME}_MergedPoint/${FILTERED_NAME}.wig
done
# wigToBigWig and fetchChrSizes are from UCSC
/usr/local/bin/wigToBigWig ${FILTERED_NAME}_MergedPoint/${FILTERED_NAME}.wig /san1_data/genomic_data/mouse/NCBI_37/$REFNAME.chrom.sizes ${FILTERED_NAME}.bw
cd $CURRENT
