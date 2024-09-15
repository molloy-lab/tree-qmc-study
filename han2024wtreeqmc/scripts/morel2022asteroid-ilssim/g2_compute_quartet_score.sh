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
MTHD=$5

cd $OUTDIR

if [ ! -e $GTRE_FILE ]; then
    echo "ERROR: $GTRE_FILE does not exist!"
    exit
fi

if [ $MTHD == "true_stree" ]; then
    INPUT=$STRE_TRUE
else
    INPUT="$MTHD.tre"
fi

QSUP="qsupp-wn"
OUTPUT="${MTHD}_${QSUP}"
if [ ! -e $OUTPUT.tre ]; then
    # Run ASTRAL
    $ASTER --scoring -u 2 -t 16 \
           -i $GTRE_FILE \
           -c $INPUT \
           -o $OUTPUT.tre \
           &> $OUTPUT.log
    MYQS="$(grep "Score:" $OUTPUT.log | awk '{print $2}')"
    echo "$MYMODL,$MTHD,$QSUP,$MYQS" > ${OUTPUT}_quartet_score.csv
fi

