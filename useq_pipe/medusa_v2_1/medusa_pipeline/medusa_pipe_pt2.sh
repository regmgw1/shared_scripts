#!/bin/bash

# part 2 of medusa pipeline - written by Gareth Wilson
# runs MEDIPS qc and generates wig tracks for total plus both strands. Wig tracks converted to bigwig. Also cleans up unnecessary files.
# requires installation of R (www.r-project.org/) and the bioconductor package MEDIPS (http://www.bioconductor.org/packages/2.7/bioc/html/MEDIPS.html) plus wigToBigWig (http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/)


NAME=`echo $1 |sed 's/\.sam$//'`
SPECIES=$2
REFNAME=$3
WINDOW_SIZE=$4
PATH2INPUT=$5
PATH2CHROM=$6
PATH2OUTPUT=$7
PATH2SCRIPTS=$8

#run the script in R using Rscript

$PATH2SCRIPTS/R_scripts/medips_qc.R $NAME $SPECIES $REFNAME $WINDOW_SIZE 0 $PATH2OUTPUT 

#SECTION TO GENERATE SINGLE STRANDED WIG FILES
grep '	+' $PATH2OUTPUT/$NAME.bed >$PATH2OUTPUT/${NAME}_for.bed
$PATH2SCRIPTS/R_scripts/medips_qc.R ${NAME}_for $SPECIES $REFNAME $WINDOW_SIZE 1 $PATH2OUTPUT
grep '	-' $PATH2OUTPUT/$NAME.bed >$PATH2OUTPUT/${NAME}_rev.bed
$PATH2SCRIPTS/R_scripts/medips_qc.R ${NAME}_rev $SPECIES $REFNAME $WINDOW_SIZE 1 $PATH2OUTPUT

# check to see if rpkm.wig file has content. if so delete filter.sam and .medips files
FILE=$PATH2OUTPUT/${NAME}_rpkm.wig
if [[ -s $FILE ]] ; then
echo "wigToBigWig -clip $PATH2OUTPUT/${NAME}_rpkm.wig $PATH2CHROM $PATH2OUTPUT/${NAME}_rpkm.bw"
wigToBigWig -clip $PATH2OUTPUT/${NAME}_rpkm.wig $PATH2CHROM $PATH2OUTPUT/${NAME}_rpkm.bw
gzip $PATH2OUTPUT/${NAME}_rpkm.wig
gzip $PATH2OUTPUT/${NAME}_for_rpkm.wig
gzip $PATH2OUTPUT/${NAME}_rev_rpkm.wig
echo "$FILE has data. Removing $PATH2INPUT/${NAME}_filter.sam, $PATH2OUTPUT/${NAME}_for.bed and $PATH2OUTPUT/${NAME}_rev.bed"
rm $PATH2INPUT/${NAME}_filter.sam
rm $PATH2INPUT/${NAME}_for.bed
rm $PATH2INPUT/${NAME}_rev.bed
else
echo "$FILE is empty. $PATH2INPUT/${NAME}_filter.sam will not be deleted"
fi ;
