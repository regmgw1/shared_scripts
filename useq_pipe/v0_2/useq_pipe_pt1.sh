#!/bin/bash

FASTQC_HOME=/usr/local/bioinf/FastQC
REFERENCE=/san1_data/genomic_data/human/GRCh37/human_GRCh37.fa
PATH2SCRIPTS=/home/POGBDOM/regmgw1/useq_pipe/v0_2
NAME=`echo $1 |sed 's/\.sam$//'`
SPECIES=$2
READ_LENGTH=$3
PATH2INPUT=$4
PATH2OUTPUT=$5

#mkdir $PATH2OUTPUT
# run fastqc on bam file
#$FASTQC_HOME/fastqc $PATH2INPUT/$NAME.sorted.bam
# convert to BED6
#$PATH2SCRIPTS/sam2bedgff.pl $PATH2INPUT/${NAME}_filter.sam $NAME $SPECIES $READ_LENGTH 0 1 $PATH2OUTPUT/$NAME.bed 
# run MEDIPs
cut -f 1,2,3,6 $PATH2OUTPUT/$NAME.bed >$PATH2OUTPUT/$NAME.medips
#run the script in R using Rscript
$PATH2SCRIPTS/medips_summary.R $NAME $SPECIES $PATH2SCRIPTS $PATH2OUTPUT 
gzip $PATH2OUTPUT/${NAME}_rpm.wig
gzip $PATH2OUTPUT/${NAME}_for_rpm.wig
gzip $PATH2OUTPUT/${NAME}_rev_rpm.wig
#rm $PATH2OUTPUT/$NAME.medips

# create sample.list file for use in pt5
cat $NAME >>$PATH2OUTPUT/sample.list
