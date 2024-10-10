#!/bin/bash

GTRE=$1

# Define directories
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
INDIR="$PROJDIR/data/cloutier2019whole/gene-trees"
OUTDIR="$PROJDIR/data/cloutier2019whole/species-trees-analysis"
ASTEROID="$GROUPDIR/group/software-compiled-on-EPYC-7313/Asteroid-nompi/Asteroid/build/bin/asteroid"

cd $OUTDIR

MTHD="asteroid"

# Asteroid
STRE="${MTHD}_${GTRE}"
if [ ! -e $STRE ]; then
    $ASTEROID -i $INDIR/$GTRE -p $STRE &> $STRE.log
fi

