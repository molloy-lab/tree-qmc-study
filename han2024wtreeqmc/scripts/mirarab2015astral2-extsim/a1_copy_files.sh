#!/bin/bash

NTAX=$1
HGHT=$2
RATE=$3
REPL=$4
NGEN=$5

MODL="model.$NTAX.$HGHT.$RATE"
MYMODL="$NTAX,$HGHT,$RATE,$REPL,NA,$NGEN"

# Define directories and files
LABDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$LABDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"

# Define software files
COMPARE="$PROJDIR/tools/compare_two_trees.py"

# Define input files
INDIR="$LABDIR/group/data/mirarab2015astral2/$MODL/$REPL"
STRE="s_tree.trees"
GTRE="estimatedgenetre"
CAML="concatenatedtree.genes${NGEN}" 

# Define output directory
OUTDIR="$PROJDIR/data/mirarab2015astral2-extsim/$MODL/$REPL"
cd $OUTDIR

cp $INDIR/$STRE $STRE
cp $INDIR/$GTRE $GTRE

if [ -e $INDIR/$CAML ]; then
    MYMTHD="caml"
    MYSTRE="$MYMTHD_${NGEN}gen"
    cp $INDIR/$CAML $MYSTRE.tre
    MYERRR=$(python3 $COMPARE -t1 $STRE -t2 $MYSTRE.tre)
    echo "$MYMODL,$MYMTHD,$MYERRR" > ${MYSTRE}_species_tree_error.csv
fi

