#!/bin/bash

# part 4 of medusa pipeline - written by Gareth Wilson
# responsible for running the dmr annotation scripts. many prerequisites for this to work correctly.
# requires bedtools to be installed and in path. Also needs all features to be written in GFF file. Contact author for example of file and the required directory structure. 

# $1 = dmr_output_folder
# $2 = path2genome_fa
# $3 = path2genes_gff
# $4 = upstream nearest gene threshold (e.g. 10000)
# $5 = downstream nearest gene threshold (e.g. 5000)
# $6 = path2features_dir
# $7 = path2features.list
# $8 = path2repeats.list
# $9 = intersect threshold (e.g. 0.1 for 10% overlap, 0.000000000000001 for single base overlap)
# $10 = path2scripts
# $11 = path2output
# $12 = pathfiltered 

GFF_OR_BED=1
CURRENT=`pwd`
cd $1
DATE=`date`
WORKING=`pwd`
mkdir ${11}
echo $WORKING
echo "$DATE
dmr_output_folder = $1
path2genome_fa = $2
path2genes_gff = $3
upstream nearest gene threshold (e.g. 10000) = $4
downstream nearest gene threshold (e.g. 5000) = $5
path2features_dir = $6
path2features.list = $7
path2repeats.list = $8
intersect threshold (e.g. 0.1 for 10% overlap, 0.000000000000001 for single base overlap) = $9
path2scripts = ${10}
path2output = ${11}
path2filtered = ${12}
path2celllist = ${13}
species = ${14}
path2chromlist = ${15}
pvalue= ${16}" >${11}/annotation.log

echo "perl ${10}/perl_scripts/peak_cg_bedtools_v2.pl $2 $WORKING $WORKING/dmrFiles.txt $3 ${15} $4 $5 ${11}"
perl ${10}/perl_scripts/peak_cg_bedtools_v2.pl $2 $WORKING $WORKING/dmrFiles.txt $3 ${15} $4 $5 ${11}

echo "${10}/R_scripts/medusaVolcano.R ${16} $WORKING ${11}"
${10}/R_scripts/medusaVolcano.R ${16} $WORKING ${11}

echo "perl ${10}/perl_scripts/peaks_in_features_bedtools_wrap_v2.pl $WORKING $WORKING/dmrFiles.txt $7 $8 $6 $9 $GFF_OR_BED ${10}/perl_scripts ${11} ${13} ${14}"
perl ${10}/perl_scripts/peaks_in_features_bedtools_wrap_v2.pl $WORKING $WORKING/dmrFiles.txt $7 $8 $6 $9 $GFF_OR_BED ${10}/perl_scripts ${11} ${13} ${14}

# uses fixed overlap threshold of 0.51
echo "perl ${10}/perl_scripts/intersect2fisher.pl $WORKING $WORKING/dmrFiles.txt $WORKING/thresholdWindows.bed $6 $7 ${10}/R_scripts ${11}"
perl ${10}/perl_scripts/intersect2fisher.pl $WORKING $WORKING/dmrFiles.txt $WORKING/thresholdWindows.bed $6 $7 ${10}/R_scripts ${11}
cd $CURRENT

