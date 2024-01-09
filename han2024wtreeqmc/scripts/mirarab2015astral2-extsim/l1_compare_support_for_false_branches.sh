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
COMPARE="$PROJDIR/tools/compare_support_for_false_branches.py"

# Define input files
DATADIR="$PROJDIR/data/mirarab2015astral2-extsim/$MODL/$REPL"

MTHD1="aster_h_t16"
MTHD2="wtreeqmc_wh_n2"

STRE_TRUE="$DATADIR/true_stree_${SUPP}_${NGEN}gen-ah-annotated.tre"
STRE_EST1="$DATADIR/${MTHD1}_${SUPP}_${NGEN}gen-ah-annotated.tre"
STRE_EST2="$DATADIR/${MTHD2}_${SUPP}_${NGEN}gen-ah-annotated.tre"
OUTF="$DATADIR/false_branches_${MTHD1}_vs_${MTHD2}_for_${SUPP}_${NGEN}gen-ah-annotated.csv"

for TREE in $STRE_TRUE $STRE_EST1 $STRE_EST2; do
    if [ ! -e $TREE ]; then
        echo "$TREE is missing!"
	exit
    fi

    if [ -z $(grep ";" $TREE) ]; then
        echo "$TREE is empty!"
	exit
    fi
done

if [ ! -e $OUTF ]; then
    python3 $COMPARE -t $STRE_TRUE -e1 $STRE_EST1 -e2 $STRE_EST2 -p $MYMODL &> $OUTF
fi

