#!/bin/bash

# part 3 of medusa pipeline - written by Gareth Wilson
# creates point directory folders for use in pt4. These directories can also be used in the Useq tools QCseqs (http://useq.sourceforge.net/cmdLnMenus.html#QCSeqs) - though this is as yet not part of medusa.

USEQ_HOME=$4
REFNAME=$3
FILTERED_NAME=`echo $1 |sed 's/\.sam$//'`
PATH2INPUT=$2

CURRENT=`pwd`
cd $PATH2INPUT
# convert to USeq binary mid-point format
# creates folder ${NAME}_Point
java -Xmx6G -jar $USEQ_HOME/Tag2Point -f $FILTERED_NAME.bed -v $REFNAME -b -i 5 
cd $CURRENT
