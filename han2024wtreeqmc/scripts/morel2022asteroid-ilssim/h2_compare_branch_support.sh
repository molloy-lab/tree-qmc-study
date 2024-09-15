#!/bin/bash

# Define software and utilities
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
COMPARE="$PROJDIR/tools/compare_branch_support.py"

MYMODL=$1
OUTDIR=$2

cd $OUTDIR

MTHD1="aster_v1.16.3.4"
MTHD2="wtreeqmc_wf_n2_refined"

QSUP="qsupp-wn"
STRE_TRUE="true_stree_${QSUP}.tre"
STRE_EST1="${MTHD1}_${QSUP}.tre"
STRE_EST2="${MTHD2}_${QSUP}.tre"
OUTF="branch_support_${MTHD1}_vs_${MTHD2}_for_${QSUP}.csv"

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

