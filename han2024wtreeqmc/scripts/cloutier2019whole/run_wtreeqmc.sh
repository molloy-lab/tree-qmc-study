#!/bin/bash

GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
DATADIR="$PROJDIR/data/cloutier2019whole/species-trees-analysis"
WTREEQMC="$PROJDIR/software/TREE-QMC/tree-qmc"

cd $DATADIR

# TREE-QMC w/ hybrid weighting scheme
MTHD="wtreeqmc_hybrid_n2"
GTRE_FILE="UCE_minus_105_abayes_gene_trees_sorted_plus_bad105.tre"
STRE_FILE="${MTHD}_${GTRE_FILE}"
if [ ! -e $STRE_FILE ]; then
    $WTREEQMC -w h -n 0.333 -x 1 -d 0.333 --root galGal -i $GTRE_FILE -o $STRE_FILE &> $STRE_FILE.log
fi

# TREE-QMC w/ support weighting scheme
MTHD="wtreeqmc_support_n2"
GTRE_FILE="UCE_minus_105_abayes_gene_trees_sorted_plus_bad105.tre"
STRE_FILE="${MTHD}_${GTRE_FILE}"
if [ ! -e $STRE_FILE ]; then
    $WTREEQMC -w s -n 0.333 -x 1 -d 0.333 -i $GTRE_FILE -o $STRE_FILE &> $STRE_FILE.log
fi

# TREE-QMC w/ length weighting scheme
MTHD="wtreeqmc_length_n2"
GTRE_FILE="UCE_minus_105_abayes_gene_trees_sorted_plus_bad105.tre"
STRE_FILE="${MTHD}_${GTRE_FILE}"
if [ ! -e $STRE_FILE ]; then
    $WTREEQMC -w l -i $GTRE_FILE -o $STRE_FILE &> $STRE_FILE.log
fi

# TREE-QMC w/ no weighting scheme
MTHD="wtreeqmc_none_n2"
GTRE_FILE="UCE_minus_105_abayes_gene_trees_sorted_plus_bad105.tre"
STRE_FILE="${MTHD}_${GTRE_FILE}"
if [ ! -e $STRE_FILE ]; then
    $WTREEQMC -w n -i $GTRE_FILE -o $STRE_FILE &> $STRE_FILE.log
fi

# TREE-QMC w/ no weighting scheme - with WT Tinamou and Chicken removed from bad 105 genes
MTHD="wtreeqmc_none_n2"
GTRE_FILE="UCE_minus_105_abayes_gene_trees_sorted_plus_bad105_without_galGal_tinGut.tre"
STRE_FILE="${MTHD}_${GTRE_FILE}"
if [ ! -e $STRE_FILE ]; then
    $WTREEQMC -w n -i $GTRE_FILE -o $STRE_FILE &> $STRE_FILE.log
fi

# TREE-QMC w/ no weighting scheme - after filtering bad 105 genes
MTHD="wtreeqmc_none_n2"
GTRE_FILE="UCE_minus_105_abayes_gene_trees_sorted.tre"
STRE_FILE="${MTHD}_${GTRE_FILE}"
if [ ! -e $STRE_FILE ]; then
    $WTREEQMC -w n -i $GTRE_FILE -o $STRE_FILE &> $STRE_FILE.log
fi

