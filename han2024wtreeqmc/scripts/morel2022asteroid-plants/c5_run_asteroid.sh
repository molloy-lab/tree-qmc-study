#!/bin/bash

# Define software and utilities
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
COMPARE="$PROJDIR/tools/compare_two_trees.py"
ASTEROID="$GROUPDIR/group/software-compiled-on-EPYC-7313/Asteroid-nompi/Asteroid/build/bin/asteroid"

MYINFO=$1
DATADIR=$2
GTRE_FILE=$3

cd $DATADIR

if [ ! -e $GTRE_FILE ]; then
    echo "ERROR: $GTRE_FILE does not exist!"
    exit
fi

MYMTHD="asteroid"
MYSTRE="$MYMTHD.bestTree.newick"
if [ -e $MYSTRE ]; then
    X=$(grep ";" $MYSTRE | wc -l)
    if [ $X -eq 0 ]; then
        rm ${MYMTHD}*
    fi
fi

if [ ! -e $MYSTRE ]; then
    MYTIME="$(time ($ASTEROID -i $GTRE_FILE \
                              -p $MYMTHD \
                              &> $MYMTHD.log) 2>&1 1>/dev/null)"

    uname -a > ${MYMTHD}_node_info.csv
    MYNODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

    MYSECS="$(echo $MYTIME | awk '{print $2","$4","$6}')"
    echo "$MYINFO,$MYMTHD,$MYNODE,$MYSECS" > ${MYMTHD}_runtime.csv
fi

