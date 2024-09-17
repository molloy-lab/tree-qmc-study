#!/bin/bash

# Define directories
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
DATADIR="$PROJDIR/data/cloutier2019whole/species-trees-analysis"
ASTEROID="$GROUPDIR/group/software-compiled-on-EPYC-7313/Asteroid-nompi/Asteroid/build/bin/asteroid"

cd $DATADIR

MTHD="asteroid"

# Asteroid
GTRE_FILE="UCE_minus_105_abayes_gene_trees_sorted_plus_bad105.tre"
STRE_FILE="${MTHD}_${GTRE_FILE}"
if [ ! -e $STRE_FILE ]; then
    $ASTEROID -i $GTRE_FILE -p $STRE_FILE &> $STRE_FILE.log
fi

# Asteroid - with WT Tinamou and Chicken removed from bad 105 genes
GTRE_FILE="UCE_minus_105_abayes_gene_trees_sorted_plus_bad105_without_galGal_tinGut.tre"
STRE_FILE="${MTHD}_${GTRE_FILE}"
if [ ! -e $STRE_FILE ]; then
    $ASTEROID -i $GTRE_FILE -p $STRE_FILE &> $STRE_FILE.log
fi

# Asteroid - after filtering bad 105 genes
GTRE_FILE="UCE_minus_105_abayes_gene_trees_sorted.tre"
STRE_FILE="${MTHD}_${GTRE_FILE}"
if [ ! -e $STRE_FILE ]; then
    $ASTEROID -i $GTRE_FILE -p $STRE_FILE &> $STRE_FILE.log
fi

