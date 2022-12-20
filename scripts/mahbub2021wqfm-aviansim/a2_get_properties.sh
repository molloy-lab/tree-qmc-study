#!/bin/bash

# Evaluate model conditions for avian simulated data set

SCAL=$1
NGEN="1000"
NBPS="500"
REPL=$2

#for REPL in `seq -f "R%g" 1 20`; do

echo "Working on $REPL..."
PREFIX="$SCAL,$NGEN,$NBPS,$REPL"

# Define directories
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJECTDIR="$GROUPDIR/ekmolloy/tree-qmc-study"
TOOLDIR="$PROJECTDIR/tools/"

# Define tools
COMPARE_LISTS="$PROJECTDIR/tools/compare_two_tree_lists.py"
COMPARE_TREE_TO_LIST="$PROJECTDIR/tools/compare_tree_to_tree_list.py"

# Define input files
OUTDIR="$PROJECTDIR/data/mahbub2021wqfm-aviansim"
STRE_TRUE="$GROUPDIR/group/data/mahbub2021wqfm/48-taxon-avian-simulated/True.tree"
GTRE_TRUE="$OUTDIR/$SCAL-$NGEN-true/$REPL/all_gt.tre"
GTRE_ESTI="all_gt.tre"

cd $OUTDIR/$SCAL-$NGEN-$NBPS/$REPL

# Difference - true species tree and true gene trees
OUTF="true_species_tree_vs_true_gene_trees.csv"
if [ ! -e $OUTF ]; then
/opt/local/stow/Python3-3.8.1/bin/python3 $COMPARE_TREE_TO_LIST \
    -t $STRE_TRUE \
    -l $GTRE_TRUE \
    -p $PREFIX \
    -o $OUTF
fi


# Difference - true species tree and estimated gene trees
OUTF="true_species_tree_vs_estimated_gene_trees.csv"
if [ ! -e $OUTF ]; then
/opt/local/stow/Python3-3.8.1/bin/python3 $COMPARE_TREE_TO_LIST \
    -t $STRE_TRUE \
    -l $GTRE_ESTI \
    -p $PREFIX \
    -o $OUTF
fi


# Difference - true and estimate gene trees
OUTF="true_vs_estimated_gene_trees.csv"
if [ ! -e $OUTF ]; then
/opt/local/stow/Python3-3.8.1/bin/python3 $COMPARE_LISTS \
   -l1 $GTRE_TRUE \
   -l2 $GTRE_ESTI \
   -p $PREFIX \
   -o $OUTF
fi

#done

