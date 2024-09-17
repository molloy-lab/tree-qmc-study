#!/bin/bash

# Define directories
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
DATADIR="$PROJDIR/data/cloutier2019whole/species-trees-analysis"
ASTERH="$GROUPDIR/group/software-compiled-on-EPYC-7313/ASTER/bin/astral-hybrid"
# --mode	Integer	1	1: hybrid weighting, 2: support only, 3: length only, 4: unweighted

cd $DATADIR

# ASTER w/ hybrid weighting scheme
MTHD="waster_hybrid"
GTRE_FILE="UCE_minus_105_abayes_gene_trees_sorted_plus_bad105.tre"
STRE_FILE="${MTHD}_${GTRE_FILE}"
if [ ! -e $STRE_FILE ]; then
    $ASTERH --mode 1 --bayes -i $GTRE_FILE -o $STRE_FILE &> $STRE_FILE.log
fi

# ASTER w/ support weighting scheme
MTHD="waster_support"
GTRE_FILE="UCE_minus_105_abayes_gene_trees_sorted_plus_bad105.tre"
STRE_FILE="${MTHD}_${GTRE_FILE}"
if [ ! -e $STRE_FILE ]; then
    $ASTERH --mode 2 --bayes -i $GTRE_FILE -o $STRE_FILE &> $STRE_FILE.log
fi

# ASTER w/ length weighting scheme
MTHD="waster_length"
GTRE_FILE="UCE_minus_105_abayes_gene_trees_sorted_plus_bad105.tre"
STRE_FILE="${MTHD}_${GTRE_FILE}"
if [ ! -e $STRE_FILE ]; then
    $ASTERH --mode 3 -i $GTRE_FILE -o $STRE_FILE &> $STRE_FILE.log
fi

# ASTER w/ no weighting scheme
MTHD="waster_none"
GTRE_FILE="UCE_minus_105_abayes_gene_trees_sorted_plus_bad105.tre"
STRE_FILE="${MTHD}_${GTRE_FILE}"
if [ ! -e $STRE_FILE ]; then
    $ASTERH --mode 4 -i $GTRE_FILE -o $STRE_FILE &> $STRE_FILE.log
fi

# ASTER w/ no weighting scheme - with WT Tinamou and Chicken removed from bad 105 genes
MTHD="waster_none"
GTRE_FILE="UCE_minus_105_abayes_gene_trees_sorted_plus_bad105_without_galGal_tinGut.tre"
STRE_FILE="${MTHD}_${GTRE_FILE}"
if [ ! -e $STRE_FILE ]; then
    $ASTERH --mode 4 -i $GTRE_FILE -o $STRE_FILE &> $STRE_FILE.log
fi

# ASTER w/ no weighting scheme - after filtering bad 105 genes
MTHD="waster_none"
GTRE_FILE="UCE_minus_105_abayes_gene_trees_sorted.tre"
STRE_FILE="${MTHD}_${GTRE_FILE}"
if [ ! -e $STRE_FILE ]; then
    $ASTERH --mode 4 -i $GTRE_FILE -o $STRE_FILE &> $STRE_FILE.log
fi

