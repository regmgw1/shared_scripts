#!/bin/bash

# part 3 of medusa pipeline - written by Gareth Wilson
# calls dmrs using medips
# requires installation of R (www.r-project.org/) and the bioconductor package MEDIPS (http://www.bioconductor.org/packages/2.7/bioc/html/MEDIPS.html)


TREAT=$1
CONTROL=$2
SPECIES=$3
REFNAME=$4
READ_DEPTH=$5
DMR_SIZE=$6
PVALUE=$7
SEX=$8
PATH2INPUT=$9
PATH2OUTPUT=${10}
PATH2SCRIPTS=${11}

mkdir $PATH2OUTPUT
#run the script in R using Rscript

echo "$PATH2SCRIPTS/R_scripts/medips_dmrs.R $TREAT $CONTROL $SPECIES $REFNAME $READ_DEPTH $DMR_SIZE $PVALUE $SEX $PATH2INPUT $PATH2OUTPUT"
$PATH2SCRIPTS/R_scripts/medips_dmrs.R $TREAT $CONTROL $SPECIES $REFNAME $READ_DEPTH $DMR_SIZE $PVALUE $SEX $PATH2INPUT $PATH2OUTPUT 


