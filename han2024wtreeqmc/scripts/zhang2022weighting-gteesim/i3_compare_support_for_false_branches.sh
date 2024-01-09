#!/bin/bash

REPL=$1
NBPS=$2
SUPP=$3
NGEN=$4

MYMODL="$REPL,$NBPS,$SUPP,$NGEN"

# Define directories
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"

# Define software
COMPARE="$PROJDIR/tools/compare_support_for_false_branches.py"

# Define input files
DATADIR="$PROJDIR/data/zhang2022weighting-gteesim/S100/$REPL"

cd $DATADIR

MTHD1="aster_h_t16"
MTHD2="wtreeqmc_wh_n2"

STRE_TRUE="$DATADIR/true_stree_${SUPP}_${NBPS}bps_${NGEN}gen-ah-annotated.tre"
STRE_EST1="$DATADIR/${MTHD1}_${SUPP}_${NBPS}bps_${NGEN}gen-ah-annotated.tre"
STRE_EST2="$DATADIR/${MTHD2}_${SUPP}_${NBPS}bps_${NGEN}gen-ah-annotated.tre"
OUTF="$DATADIR/false_branches_${MTHD1}_vs_${MTHD2}_for_${SUPP}_${NBPS}bps_${NGEN}gen-ah-annotated.csv"

for TREE in $STRE_TRUE $STRE_EST1 $STRE_EST2; do
    if [ ! -e $TREE ]; then
        echo "$DATADIR/$TREE is missing!"
	exit
    fi

    if [ -z $(grep ";" $TREE) ]; then
        echo "$DATADIR/$TREE is empty!"
	exit
    fi
done

if [ ! -e $OUTF ]; then
    python3 $COMPARE -t $STRE_TRUE -e1 $STRE_EST1 -e2 $STRE_EST2 &> $OUTF
fi

