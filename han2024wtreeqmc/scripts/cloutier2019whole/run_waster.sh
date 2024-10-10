#!/bin/bash

GTRE=$1

# Define directories
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
INDIR="$PROJDIR/data/cloutier2019whole/gene-trees"
OUTDIR="$PROJDIR/data/cloutier2019whole/species-trees-analysis"
ASTERH="$GROUPDIR/group/software-compiled-on-EPYC-7313/ASTER/bin/astral-hybrid"
# --mode	Integer	1	1: hybrid weighting, 2: support only, 3: length only, 4: unweighted

cd $OUTDIR
pwd

# ASTER w/ hybrid weighting scheme
MTHD="waster_hybrid"
STRE="${MTHD}_${GTRE}"
if [ ! -e $STRE ]; then
    $ASTERH --mode 1 --bayes -i $INDIR/$GTRE -o $STRE &> $STRE.log
fi

# ASTER w/ support weighting scheme
MTHD="waster_support"
STRE="${MTHD}_${GTRE}"
if [ ! -e $STRE ]; then
    $ASTERH --mode 2 --bayes -i $INDIR/$GTRE -o $STRE &> $STRE.log
fi

# ASTER w/ length weighting scheme
MTHD="waster_length"
STRE="${MTHD}_${GTRE}"
if [ ! -e $STRE ]; then
    $ASTERH --mode 3 -i $INDIR/$GTRE -o $STRE &> $STRE.log
fi

# ASTER w/ no weighting scheme
MTHD="waster_none"
STRE="${MTHD}_${GTRE}"
if [ ! -e $STRE ]; then
    $ASTERH --mode 4 -i $INDIR/$GTRE -o $STRE &> $STRE.log
fi

