#!/bin/bash

GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJECTDIR="$GROUPDIR/ekmolloy/species-to-tumor-study"

BUILDTREE="$PROJECTDIR/tools/construct_tree_from_conflict_free_matrix.py"
COMPARE="$PROJECTDIR/tools/compare_two_trees.py"
CALCERROR="$PROJECTDIR/tools/calculate_mut_pair_error.py"
CALCWHUNTRESS="$PROJECTDIR/tools/calculate_mut_pair_error_with_huntress.py"

MYMODL=$1
DATDIR=$2

cd $DATDIR

pwd

MYMTHD="huntress_v0.1.2.0_default"

TRUE_MUTS="ground.CFMatrix"
ESTI_MUTS="$MYMTHD.CFMatrix"
FOUT_MUTS="${MYMTHD}_mut_pair_error.csv"
if [ -e $ESTI_MUTS ] && [ ! -e $FOUT_MUTS ]; then    
    MYERRR=$(/opt/local/stow/Python3-3.8.1/bin/python3 $CALCERROR -t $TRUE_MUTS -e $ESTI_MUTS)
    echo "$MYMODL,$MYMTHD,$MYERRR" > $FOUT_MUTS
fi

TRUE_TREE="ground.cell_lineage_tree.nwk"
ESTI_TREE="$MYMTHD.cell_lineage_tree.nwk"
FOUT_TREE="${MYMTHD}_cell_lineage_tree_error.csv"

if [ -e $ESTI_MUTS ] && [ ! -e $ESTI_TREE ]; then
        /opt/local/stow/Python3-3.8.1/bin/python3 $BUILDTREE \
            -i $ESTI_MUTS \
            -o $ESTI_TREE
fi

if [ -e $ESTI_TREE ] && [ ! -e $FOUT_TREE ]; then
    MYERRR=$(/opt/local/stow/Python3-3.8.1/bin/python3 $COMPARE -t1 $TRUE_TREE -t2 $ESTI_TREE)
    echo "$MYMODL,$MYMTHD,$MYERRR" > $FOUT_TREE
fi

