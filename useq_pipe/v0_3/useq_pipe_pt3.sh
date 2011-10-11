#!/bin/bash

USEQ_HOME=/usr/local/bioinf/USeq_6.8/Apps
PERL_SCRIPTS=$5
REFNAME=$4
PATH2GENOME=$3
FILTERED_NAME=`echo $1 |sed 's/\.sam$//'`
PATH2INPUT=$2

CURRENT=`pwd`
cd $PATH2INPUT
# convert to USeq binary mid-point format
# creates folder ${NAME}_Point
java -Xmx6G -jar $USEQ_HOME/Tag2Point -f $FILTERED_NAME.bed -v $REFNAME -b -i 5 
cd $CURRENT
