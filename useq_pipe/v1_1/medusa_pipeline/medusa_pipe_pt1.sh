#!/bin/bash

# part 1 of medusa pipeline - written by Gareth Wilson
# performs final filtering step and some medip specific QC
# requires installation of fastqc (www.bioinformatics.bbsrc.ac.uk/projects/fastqc/), R (www.r-project.org/) and the bioconductor package MEDIPS (http://www.bioconductor.org/packages/2.7/bioc/html/MEDIPS.html)


NAME=`echo $1 |sed 's/\.sam$//'`
SPECIES=$2
READ_LENGTH=$3
PATH2INPUT=$4
PATH2OUTPUT=$5
PATH2SCRIPTS=$6

# run fastqc on bam file
fastqc $PATH2INPUT/$NAME.sorted.bam
# convert to BED
echo "$PATH2SCRIPTS/perl_scripts/sam2bedgff.pl $PATH2INPUT/${NAME}_filter.sam $NAME $SPECIES $READ_LENGTH 0 1 $PATH2OUTPUT/$NAME.bed";
$PATH2SCRIPTS/perl_scripts/sam2bedgff.pl $PATH2INPUT/${NAME}_filter.sam $NAME $SPECIES $READ_LENGTH 0 1 $PATH2OUTPUT/$NAME.bed 
# generate input file for MEDIPS
cut -f 1,2,3,6 $PATH2OUTPUT/$NAME.bed >$PATH2OUTPUT/$NAME.medips
#run the script in R using Rscript
$PATH2SCRIPTS/R_scripts/medips_summary.R $NAME $SPECIES $PATH2SCRIPTS $PATH2OUTPUT 

# check to see if rpm.wig file has content. if so delete filter.sam and .medips files
FILE=$PATH2OUTPUT/${NAME}_rpm.wig
if [[ -s $FILE ]] ; then
wigToBigWig $PATH2OUTPUT/${NAME}_rpm.wig chrom.sizes $PATH2OUTPUT/${NAME}_rpm.bw
gzip $PATH2OUTPUT/${NAME}_rpm.wig
gzip $PATH2OUTPUT/${NAME}_for_rpm.wig
gzip $PATH2OUTPUT/${NAME}_rev_rpm.wig
echo "$FILE has data. Removing $PATH2INPUT/${NAME}_filter.sam and $PATH2OUTPUT/$NAME.medips"
rm $PATH2INPUT/${NAME}_filter.sam
rm $PATH2OUTPUT/$NAME.medips
else
echo "$FILE is empty. $PATH2INPUT/${NAME}_filter.sam and $PATH2OUTPUT/$NAME.medips will not be deleted"
fi ;
