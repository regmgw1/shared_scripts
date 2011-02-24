#!/bin/bash

# $1 = useq_output_folder
# $2 = path2genome_dir
# $3 = path2genes_gff
# $4 = upstream nearest gene threshold (e.g. 10000)
# $5 = downstream nearest gene threshold (e.g. 5000)
# $6 = path2features_dir
# $7 = path2features.list
# $8 = path2repeats.list
# $9 = intersect threshold (e.g. 0.1 for 10% overlap, 0.000000000000001 for single base overlap)
# $10 = path2scripts
# $11 = path2output 

GFF_OR_BED=0
CURRENT=`pwd`
cd $1
ls -p|grep /|grep FDR >peakDir.list
DATE=`date`
WORKING=`pwd`
mkdir ${11}
echo $WORKING
echo "$DATE
useq_output_folder = $1
path2genome_dir = $2
path2genes_gff = $3
upstream nearest gene threshold (e.g. 10000) = $4
downstream nearest gene threshold (e.g. 5000) = $5
path2features_dir = $6
path2features.list = $7
path2repeats.list = $8
intersect threshold (e.g. 0.1 for 10% overlap, 0.000000000000001 for single base overlap) = $9
path2scripts = ${10}
path2output = ${11}" >${11}/annotation.log

echo "perl ${10}/peak_cg_bedtools.pl $2 $WORKING $WORKING/peakDir.list $3 $4 $5 ${11}"
perl ${10}/peak_cg_bedtools.pl $2 $WORKING $WORKING/peakDir.list $3 $4 $5 ${11}

echo "perl ${10}/peak_explorer_bedtools.pl ../sample.list ../ $WORKING $WORKING/peakDir.list ${11}"
perl ${10}/peak_explorer_bedtools.pl ../sample.list ../ $WORKING $WORKING/peakDir.list ${11}

echo "perl ${10}/peaks_in_features_bedtools_wrap.pl $WORKING $WORKING/peakDir.list $7 $8 $6 $9 $GFF_OR_BED ${11}"
perl ${10}/peaks_in_features_bedtools_wrap.pl $WORKING $WORKING/peakDir.list $7 $8 $6 $9 $GFF_OR_BED ${10} ${11}
cd $CURRENT
