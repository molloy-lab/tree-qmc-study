#!/bin/bash

# Define software and utilities
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
COMPARE="$PROJDIR/tools/compare_two_trees.py"
ASTER="$GROUPDIR/group/software-compiled-on-EPYC-7313/ASTER/bin/astral"

MYMODL=$1
OUTDIR=$2
GTRE_FILE=$3
STRE_TRUE=$4

cd $OUTDIR

if [ ! -e $GTRE_FILE ]; then
    echo "ERROR: $GTRE_FILE does not exist!"
    exit
fi

MYMTHD="aster_v1.16.3.4"
if [ -e $MYMTHD.tre ]; then
    X=$(grep ";" $MYMTHD.tre | wc -l)
    if [ $X -eq 0 ]; then
        rm ${MYMTHD}*
    fi
fi

if [ ! -e $MYMTHD.tre ]; then
    MYTIME="$(time ($ASTER -u 0 -i $GTRE_FILE \
                           -o $MYMTHD.tre \
                           &> $MYMTHD.log) 2>&1 1>/dev/null)"

    uname -a > ${MYMTHD}_node_info.csv
    MYNODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

    MYSECS="$(echo $MYTIME | awk '{print $2","$4","$6}')"
    echo "$MYMODL,$MYMTHD,$MYNODE,$MYSECS" > ${MYMTHD}_runtime.csv

    MYERRR=$(python3 $COMPARE -t1 $STRE_TRUE -t2 $MYMTHD.tre)
    echo "$MYMODL,$MYMTHD,$MYERRR" > ${MYMTHD}_species_tree_error.csv
fi

