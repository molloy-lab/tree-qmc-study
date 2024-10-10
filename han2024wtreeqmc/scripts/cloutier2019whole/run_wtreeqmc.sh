#!/bin/bash

GTRE=$1

GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
INDIR="$PROJDIR/data/cloutier2019whole/gene-trees"
OUTDIR="$PROJDIR/data/cloutier2019whole/species-trees-analysis"
WTREEQMC="$PROJDIR/software/TREE-QMC/tree-qmc"

cd $OUTDIR

# TREE-QMC w/ hybrid weighting scheme
MTHD="wtreeqmc_hybrid_n2"
STRE="${MTHD}_${GTRE}"
if [ ! -e $STRE ]; then
    $WTREEQMC -w h -n 0.333 -x 1 -d 0.333 --root galGal -i $INDIR/$GTRE -o $STRE &> $STRE.log
fi

# TREE-QMC w/ support weighting scheme
MTHD="wtreeqmc_support_n2"
STRE="${MTHD}_${GTRE}"
if [ ! -e $STRE ]; then
    $WTREEQMC -w s -n 0.333 -x 1 -d 0.333 -i $INDIR/$GTRE -o $STRE &> $STRE.log
fi

# TREE-QMC w/ length weighting scheme
MTHD="wtreeqmc_length_n2"
STRE="${MTHD}_${GTRE}"
if [ ! -e $STRE ]; then
    $WTREEQMC -w l -i $INDIR/$GTRE -o $STRE &> $STRE.log
fi

# TREE-QMC w/ no weighting scheme
MTHD="wtreeqmc_none_n2"
STRE="${MTHD}_${GTRE}"
if [ ! -e $STRE ]; then
    $WTREEQMC -w n -i $INDIR/$GTRE -o $STRE &> $STRE.log
fi

