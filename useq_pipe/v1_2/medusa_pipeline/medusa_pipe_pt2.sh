#!/bin/bash

# part 1 of medusa pipeline - written by Gareth Wilson
# performs final filtering step and some medip specific QC
# requires installation of fastqc (www.bioinformatics.bbsrc.ac.uk/projects/fastqc/), R (www.r-project.org/) and the bioconductor package MEDIPS (http://www.bioconductor.org/packages/2.7/bioc/html/MEDIPS.html)


NAME=`echo $1 |sed 's/\.sam$//'`
SPECIES=$2
WINDOW_SIZE=$3
PATH2INPUT=$4
PATH2CHROM=$5
PATH2OUTPUT=$6
PATH2SCRIPTS=$7

#grep -v '^GL\|^M\|^chrGL\|^chrM' $PATH2OUTPUT/$NAME.bed|cut -f 1,2,3,6 >$PATH2OUTPUT/$NAME.medips
#run the script in R using Rscript

$PATH2SCRIPTS/R_scripts/medips_qc.R $NAME $SPECIES $WINDOW_SIZE 0 $PATH2OUTPUT 

#SECTION TO GENERATE SINGLE STRANDED WIG FILES
grep '	+' $PATH2OUTPUT/$NAME.bed >$PATH2OUTPUT/${NAME}_for.bed
$PATH2SCRIPTS/R_scripts/medips_qc.R ${NAME}_for $SPECIES $WINDOW_SIZE 1 $PATH2OUTPUT
grep '	-' $PATH2OUTPUT/$NAME.bed >$PATH2OUTPUT/${NAME}_rev.bed
$PATH2SCRIPTS/R_scripts/medips_qc.R ${NAME}_rev $SPECIES $WINDOW_SIZE 1 $PATH2OUTPUT

# check to see if rpm.wig file has content. if so delete filter.sam and .medips files
FILE=$PATH2OUTPUT/${NAME}_rpm.wig
if [[ -s $FILE ]] ; then
echo "wigToBigWig -clip $PATH2OUTPUT/${NAME}_rpm.wig $PATH2CHROM $PATH2OUTPUT/${NAME}_rpm.bw"
wigToBigWig -clip $PATH2OUTPUT/${NAME}_rpm.wig $PATH2CHROM $PATH2OUTPUT/${NAME}_rpm.bw
gzip $PATH2OUTPUT/${NAME}_rpm.wig
gzip $PATH2OUTPUT/${NAME}_for_rpm.wig
gzip $PATH2OUTPUT/${NAME}_rev_rpm.wig
echo "$FILE has data. Removing $PATH2INPUT/${NAME}_filter.sam, $PATH2OUTPUT/${NAME}_for.bed and $PATH2OUTPUT/${NAME}_rev.bed"
rm $PATH2INPUT/${NAME}_filter.sam
rm $PATH2INPUT/${NAME}_for.bed
rm $PATH2INPUT/${NAME}_rev.bed
else
echo "$FILE is empty. $PATH2INPUT/${NAME}_filter.sam will not be deleted"
fi ;
