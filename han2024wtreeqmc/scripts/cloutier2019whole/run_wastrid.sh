#!/bin/bash

GTRE=$1

GROUPDIR="/fs/cbcb-lab/ekmolloy"
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
INDIR="$PROJDIR/data/cloutier2019whole/gene-trees"
OUTDIR="$PROJDIR/data/cloutier2019whole/species-trees-analysis"
WASTRID="$GROUPDIR/group/software/wastrid"

cd $OUTDIR

# ASTRID w/ support weighting scheme
MTHD="wastrid_support"
STRE="${MTHD}_${GTRE}"
if [ ! -e $STRE ]; then
    $WASTRID -m support -b 0.333-1 -i $INDIR/$GTRE -o $STRE &> $STRE.log
fi

# ASTRID w/ length weighting scheme
MTHD="wastrid_length"
STRE="${MTHD}_${GTRE}"
if [ ! -e $STRE ]; then
    $WASTRID -m n-length -i $INDIR/$GTRE -o $STRE &> $STRE.log
fi

# ASTRID w/ no weighting scheme
MTHD="wastrid_none"
STRE="${MTHD}_${GTRE}"
if [ ! -e $STRE ]; then
    $WASTRID -m internode -i $INDIR/$GTRE -o $STRE &> $STRE.log
fi

