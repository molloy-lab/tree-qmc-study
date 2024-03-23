#!/bin/bash

REPL=$1
NBPS=$2
SUPP=$3
NGEN=$4
MTHD=$5

MYMODL="$REPL,$NBPS,$SUPP,$NGEN,$MTHD"

# Define directories
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"

# Define software
COMPARE="$PROJDIR/tools/compare_two_trees.py"
ASTER="$GROUPDIR/group/software-compiled-on-EPYC-7313/ASTER/bin/astral"
ASTERH="$GROUPDIR/group/software-compiled-on-EPYC-7313/ASTER/bin/astral-hybrid"

# Define input files
INDIR="$GROUPDIR/group/data/zhang2022weighting-dryad/S100/$REPL"
OUTDIR="$PROJDIR/data/zhang2022weighting-gteesim/S100/$REPL"
STRE_TRUE="$INDIR/s_tree.trees"
GTRE="$INDIR/bestMLestimatedgenetree/estimatedgenetre_${NBPS}.gtr.rerooted.final.contracted.non"
GTRE_FILE="estimatedgenetre.${NBPS}.${NGEN}"

# Options
ROPTS="--bootstrap"
if [ $SUPP == "abayes" ]; then
    ROPTS="--bayes"
    GTRE="$GTRE.abayes"
fi

# Do work
cd $OUTDIR

GTRE_FILE="${GTRE}.${NGEN}"
if [ ! -e $GTRE_FILE ]; then
    head -n${NGEN} $GTRE > $GTRE_FILE
fi

if [ $MTHD == "true_stree" ]; then
    NAME="${MTHD}_${SUPP}_${NBPS}bps_${NGEN}gen"
    INPUT=$STRE_TRUE
else
    NAME="${MTHD}_${SUPP}_${NBPS}bps_${NGEN}gen"
    INPUT="$NAME.tre"
fi

if [ ! -e $INPUT ]; then
    echo "$INPUT is missing!"
    exit
fi    

if [ -z $(grep ";" $INPUT) ]; then
    echo "$INPUT is empty!"
    exit
fi

QSUP="qsupp-wh"
OUTPUT="${NAME}_${QSUP}"
if [ ! -e $OUTPUT.tre ]; then
    # Run weighted ASTRAL
    $ASTERH $ROPTS \
            --scoring -u 2 -t 16 \
            -i $GTRE_FILE \
            -c $INPUT \
            -o $OUTPUT.tre \
            &> $OUTPUT.log
    MYQS="$(grep "Score:" $OUTPUT.log | awk '{print $2}')"
    echo "$MYMODL,$QSUP,$MYQS" > ${OUTPUT}_quartet_score.csv

fi

QSUP="qsupp-wn"
OUTPUT="${NAME}_${QSUP}"
if [ ! -e $OUTPUT.tre ]; then
    # Run ASTRAL
    $ASTER --scoring -u 2 -t 16 \
            -i $GTRE_FILE \
            -c $INPUT \
            -o $OUTPUT.tre \
            &> $OUTPUT.log
    MYQS="$(grep "Score:" $OUTPUT.log | awk '{print $2}')"
    echo "$MYMODL,$QSUP,$MYQS" > ${OUTPUT}_quartet_score.csv
fi

