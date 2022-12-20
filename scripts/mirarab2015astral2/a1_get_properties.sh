#!/bin/bash

# Evaluate model conditions for ASTRAL-II paper, again...

NTAX=$1
PSIZE=$2
SRATE=$3

for REPL in `seq -f "%02g" 1 50`; do

echo "Working on $REPL..."
PREFIX="$NTAX,$PSIZE,$SRATE,$REPL"

# Define directories
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJECTDIR="$GROUPDIR/ekmolloy/tree-qmc-study"
TOOLDIR="$PROJECTDIR/tools/"

MODL="model.$NTAX.$PSIZE.$SRATE"
INDIR="$GROUPDIR/group/data/mirarab2015astral2/$MODL/$REPL"
OUTDIR="$PROJECTDIR/data/mirarab2015astral2/$MODL/$REPL"

# Define tools
COMPARE_LISTS="$PROJECTDIR/tools/compare_two_tree_lists.py"
COMPARE_TREE_TO_LIST="$PROJECTDIR/tools/compare_tree_to_tree_list.py"

# Define input files
STRE_TRUE="$INDIR/s_tree.trees"
GTRE_TRUE="$INDIR/truegenetrees"
GTRE_ESTI="$INDIR/estimatedgenetre"

# Difference - true species tree and true gene trees
OUTF="$OUTDIR/true_species_tree_vs_true_gene_trees.csv"
if [ ! -e $OUTF ]; then
/opt/local/stow/Python3-3.8.1/bin/python3 $COMPARE_TREE_TO_LIST \
    -t $STRE_TRUE \
    -l $GTRE_TRUE \
    -p $PREFIX \
    -o $OUTF
fi


# Difference - true species tree and estimated gene trees
OUTF="$OUTDIR/true_species_tree_vs_estimated_gene_trees.csv"
if [ ! -e $OUTF ]; then
/opt/local/stow/Python3-3.8.1/bin/python3 $COMPARE_TREE_TO_LIST \
    -t $STRE_TRUE \
    -l $GTRE_ESTI \
    -p $PREFIX \
    -o $OUTF
fi


# Difference - true and estimate gene trees
OUTF="$OUTDIR/true_vs_estimated_gene_trees.csv"
if [ ! -e $OUTF ]; then
/opt/local/stow/Python3-3.8.1/bin/python3 $COMPARE_LISTS \
   -l1 $GTRE_TRUE \
   -l2 $GTRE_ESTI \
   -p $PREFIX \
   -o $OUTF
fi

done
