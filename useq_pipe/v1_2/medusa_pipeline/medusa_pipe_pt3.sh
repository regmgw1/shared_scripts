#!/bin/bash

# part 1 of medusa pipeline - written by Gareth Wilson
# performs final filtering step and some medip specific QC
# requires installation of fastqc (www.bioinformatics.bbsrc.ac.uk/projects/fastqc/), R (www.r-project.org/) and the bioconductor package MEDIPS (http://www.bioconductor.org/packages/2.7/bioc/html/MEDIPS.html)


TREAT=$1
CONTROL=$2
SPECIES=$3
READ_DEPTH=$4
DMR_SIZE=$5
PVALUE=$6
PATH2INPUT=$7
PATH2OUTPUT=$8
PATH2SCRIPTS=$9

mkdir $PATH2OUTPUT
#run the script in R using Rscript

echo "$PATH2SCRIPTS/R_scripts/medips_dmrs.R $TREAT $CONTROL $SPECIES $READ_DEPTH $DMR_SIZE $PVALUE $PATH2INPUT $PATH2OUTPUT"
$PATH2SCRIPTS/R_scripts/medips_dmrs.R $TREAT $CONTROL $SPECIES $READ_DEPTH $DMR_SIZE $PVALUE $PATH2INPUT $PATH2OUTPUT 


