#!/bin/bash

NTAX=$1
HGHT=$2
RATE=$3
REPL=$4
SUPP=$5
NGEN=$6

MODL="model.$NTAX.$HGHT.$RATE"
MYMODL="$NTAX,$HGHT,$RATE,$REPL,$SUPP,$NGEN"

# Define directories
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"

# Define software
COMPARE="$PROJDIR/tools/compare_two_trees.py"
ASTERH="$GROUPDIR/group/software-compiled-on-EPYC-7313/ASTER/bin/astral-hybrid"


# Define input files
DATADIR="$PROJDIR/data/mirarab2015astral2-extsim/$MODL/$REPL"
STRE_TRUE="s_tree.trees"
GTRE="estimatedgenetre"
OPTS="-n 0.0 -x 1 -d 0.0 -u 0"
if [ $SUPP == "abayes" ]; then
    OPTS="-n 0.333 -x 1 -d 0.0 -u 0"
    GTRE="$GTRE.abayes.fixed"
fi

# Do work
cd $DATADIR

GTRE_FILE="${GTRE}.${NGEN}"
if [ ! -e $GTRE_FILE ]; then
    head -n${NGEN} $GTRE > $GTRE_FILE
fi


MYMTHD="aster_h"
MYSTRE="${MYMTHD}_${SUPP}_${NGEN}gen"
if [ ! -e $MYSTRE.tre ]; then
    MYTIME="$(time ($ASTERH $OPTS \
                         -i $GTRE_FILE \
                         -o $MYSTRE.tre \
                         &> $MYSTRE.log) 2>&1 1>/dev/null)"

    uname -a > ${MYSTRE}_node_info.csv
    MYNODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

    MYSECS="$(echo $MYTIME | awk '{print $2","$4","$6}')"
    echo "$MYMODL,$MYMTHD,$MYNODE,$MYSECS" > ${MYSTRE}_runtime.csv

    MYERRR=$(python3 $COMPARE -t1 $STRE_TRUE -t2 $MYSTRE.tre)
    echo "$MYMODL,$MYMTHD,$MYERRR" > ${MYSTRE}_species_tree_error.csv
fi
