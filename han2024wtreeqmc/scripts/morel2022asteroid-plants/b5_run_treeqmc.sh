#!/bin/bash

# Define software and utilities
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
COMPARE="$PROJDIR/tools/compare_two_trees.py"
WTREEQMC="$PROJDIR/software/TREE-QMC/tree-qmc"

MYINFO=$1
DATADIR=$2
GTRE_FILE=$3
STRE_FILE=$4

cd $DATADIR

if [ ! -e $GTRE_FILE ]; then
    echo "ERROR: $GTRE_FILE does not exist!"
    exit
fi


MYMTHD="treeqmc_wf_n2"
if [ ! -e $MYMTHD.tre ]; then
    MYTIME="$(time ($WTREEQMC -w f \
	                      -i $GTRE_FILE \
                              -o $MYMTHD.tre \
                              &> $MYMTHD.log) 2>&1 1>/dev/null)"

    uname -a > ${MYMTHD}_node_info.csv
    MYNODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

    MYSECS="$(echo $MYTIME | awk '{print $2","$4","$6}')"
    echo "$MYINFO,$MYMTHD,$MYNODE,$MYSECS" > ${MYMTHD}_runtime.csv

    MYERRR=$(/opt/local/stow/Python3-3.8.1/bin/python3 $COMPARE -t1 $STRE_FILE -t2 $MYMTHD.tre)
    echo "$MYINFO,$MYMTHD,$MYERRR" > ${MYMTHD}_species_tree_error.csv
fi

