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
WTREEQMC="$PROJDIR/software/weighted-TREE-QMC/wTREE-QMC"

# Define input files
DATADIR="$PROJDIR/data/mirarab2015astral2-extsim/$MODL/$REPL"
STRE_TRUE="s_tree.trees"
GTRE="estimatedgenetre"
OPTS="-r 0 1"
if [ $SUPP == "abayes" ]; then
    OPTS="-r 0.333 1"
    GTRE="$GTRE.abayes.fixed"
fi

# Do work
cd $DATADIR

GTRE_FILE="${GTRE}.${NGEN}"
if [ ! -e $GTRE_FILE ]; then
    head -n${NGEN} $GTRE > $GTRE_FILE
fi

MYMTHDS=( "wtreeqmc_ws_n2" \
	  "wtreeqmc_wh_n1" \
	  "wtreeqmc_wh_n0" )

for MYMTHD in ${MYMTHDS[@]}; do
    if [ $MYMTHD == "wtreeqmc_ws_n2" ]; then
        OPTS="-w s -n 2 $OPTS"
    elif [ $MYMTHD == "wtreeqmc_wh_n2" ]; then
        OPTS="-w h -n 2 $OPTS"
    elif [ $MYMTHD == "wtreeqmc_wh_n1" ]; then
        OPTS="-w h -n 1 $OPTS"
    elif [ $MYMTHD == "wtreeqmc_wh_n0" ]; then
        OPTS="-w h -n 0 $OPTS"
    else
        echo "Do not recognize $MYMTHD"
	exit
    fi

    MYSTRE="${MYMTHD}_${SUPP}_${NGEN}gen"
    if [ ! -e $MYSTRE.tre ]; then
        MYTIME="$(time ($WTREEQMC $OPTS \
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
done

