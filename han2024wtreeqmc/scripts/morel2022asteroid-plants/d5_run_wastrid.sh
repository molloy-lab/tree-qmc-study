#!/bin/bash

# Define software and utilities
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
COMPARE="$PROJDIR/tools/compare_two_trees.py"
WASTRID="$GROUPDIR/group/software/wastrid"

MYINFO=$1
DATADIR=$2
GTRE_FILE=$3

cd $DATADIR

if [ ! -e $GTRE_FILE ]; then
    echo "ERROR: $GTRE_FILE does not exist!"
    exit
fi


MYMTHD="wastrid_vanilla"
if [ ! -e $MYMTHD.tre ]; then
    MYTIME="$(time ($WASTRID --preset vanilla \
	                     -i $GTRE_FILE \
                             -o $MYMTHD.tre \
                             &> $MYMTHD.log) 2>&1 1>/dev/null)"

    uname -a > ${MYMTHD}_node_info.csv
    MYNODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

    MYSECS="$(echo $MYTIME | awk '{print $2","$4","$6}')"
    echo "$MYINFO,$MYMTHD,$MYNODE,$MYSECS" > ${MYMTHD}_runtime.csv
fi

