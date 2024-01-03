#!/bin/bash

GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
COMPARE_TREE_TO_LIST="$PROJDIR/tools/compare_tree_to_tree_list.py"
COMPARE_LISTS="$PROJDIR/tools/compare_two_tree_lists.py "

MYMODL=$1
OUTDIR=$2
STRE_TRUE=$3
GTRE_TRUE=$4
GTRE_ESTI=$5

# Diff - true species tree and true gene trees
OUTF="$OUTDIR/true_species_tree_vs_true_gene_trees.csv"
if [ ! -e $OUTF ]; then
python $COMPARE_TREE_TO_LIST \
    -t $STRE_TRUE \
    -l $GTRE_TRUE \
    -p $MYMODL &> $OUTF
fi


# Diff - true species tree and estimated gene trees
OUTF="$OUTDIR/true_species_tree_vs_estimated_gene_trees.csv"
if [ ! -e $OUTF ]; then
python $COMPARE_TREE_TO_LIST \
    -t $STRE_TRUE \
    -l $GTRE_ESTI \
    -p $MYMODL &> $OUTF
fi


# Diff - true and estimate gene trees
OUTF="$OUTDIR/true_vs_estimated_gene_trees.csv"
if [ ! -e $OUTF ]; then
python $COMPARE_LISTS \
   -l1 $GTRE_TRUE \
   -l2 $GTRE_ESTI \
   -p $MYMODL &> $OUTF
fi

