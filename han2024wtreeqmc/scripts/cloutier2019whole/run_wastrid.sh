#!/bin/bash

GROUPDIR="/fs/cbcb-lab/ekmolloy"
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
DATADIR="$PROJDIR/data/cloutier2019whole/species-trees-analysis"
WASTRID="$GROUPDIR/group/software/wastrid"

cd $DATADIR

# ASTRID w/ support weighting scheme
MTHD="wastrid_support"
GTRE_FILE="UCE_minus_105_abayes_gene_trees_sorted_plus_bad105.tre"
STRE_FILE="${MTHD}_${GTRE_FILE}"
if [ ! -e $STRE_FILE ]; then
    $WASTRID -m support -b 0.333-1 -i $GTRE_FILE -o $STRE_FILE &> $STRE_FILE.log
fi

# ASTRID w/ length weighting scheme
MTHD="wastrid_length"
GTRE_FILE="UCE_minus_105_abayes_gene_trees_sorted_plus_bad105.tre"
STRE_FILE="${MTHD}_${GTRE_FILE}"
if [ ! -e $STRE_FILE ]; then
    $WASTRID -m n-length -i $GTRE_FILE -o $STRE_FILE &> $STRE_FILE.log
fi

# ASTRID w/ no weighting scheme
MTHD="wastrid_none"
GTRE_FILE="UCE_minus_105_abayes_gene_trees_sorted_plus_bad105.tre"
STRE_FILE="${MTHD}_${GTRE_FILE}"
if [ ! -e $STRE_FILE ]; then
    $WASTRID -m internode -i $GTRE_FILE -o $STRE_FILE &> $STRE_FILE.log
fi

# ASTRID w/ no weighting scheme - with WT Tinamou and Chicken removed from bad 105 genes
MTHD="wastrid_none"
GTRE_FILE="UCE_minus_105_abayes_gene_trees_sorted_plus_bad105_without_galGal_tinGut.tre"
STRE_FILE="${MTHD}_${GTRE_FILE}"
if [ ! -e $STRE_FILE ]; then
    $WASTRID -m internode -i $GTRE_FILE -o $STRE_FILE &> $STRE_FILE.log
fi

# ASTRID w/ no weighting scheme - after filtering bad 105 genes
MTHD="wastrid_none"
GTRE_FILE="UCE_minus_105_abayes_gene_trees_sorted.tre"
STRE_FILE="${MTHD}_${GTRE_FILE}"
if [ ! -e $STRE_FILE ]; then
    $WASTRID -m internode -i $GTRE_FILE -o $STRE_FILE &> $STRE_FILE.log
fi



