#!/bin/bash

# Define software and utilities
GROUPDIR="/fs/cbcb-lab/ekmolloy"

PROJECTDIR="$GROUPDIR/ekmolloy/species-to-tumor-study"
ASTER="$GROUPDIR/group/software/ASTER/bin/astral"
COMPARE="$PROJECTDIR/tools/compare_two_trees.py"

MYMODL=$1
DATDIR=$2

#MYMTHD="aster_v1.10.2.1"
MYMTHD="aster_v1.10.2.1_wroot"

#INPUT="$DATDIR/noisy.nwk"
INPUT="$DATDIR/noisy_wroot.nwk"
CTRE_TRUE="ground.cell_lineage_tree.nwk"
CTRE_ESTI="$MYMTHD.cell_lineage_tree.nwk"
NTHREAD=16

cd $DATDIR
pwd

if [ -e $CTRE_ESTI ]; then
    X=$(grep ";" $CTRE_ESTI | wc -l)
    if [ $X -eq 0 ]; then
        rm ${MYMTHD}*
    fi
fi

if [ ! -e $CTRE_ESTI ]; then
    MYTIME="$(time ($ASTER -i $INPUT \
                           -o $CTRE_ESTI \
                           --thread $NTHREAD \
                           &> $MYMTHD.log) 2>&1 1>/dev/null)"

    uname -a > ${MYMTHD}_node_info.csv
    MYNODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

    MYSECS="$(echo $MYTIME | awk '{print $2","$4","$6}')"
    echo "$MYMODL,$MYMTHD,$MYNODE,$NTHREAD,$MYSECS" > ${MYMTHD}_runtime.csv

    MYERRR=$(/opt/local/stow/Python3-3.8.1/bin/python3 $COMPARE -t1 $CTRE_TRUE -t2 $CTRE_ESTI)
    echo "$MYMODL,$MYMTHD,$MYERRR" > ${MYMTHD}_cell_lineage_tree_error.csv
fi
