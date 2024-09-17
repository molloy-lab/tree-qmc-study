#!/bin/bash

GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
DATADIR="$PROJDIR/data/cloutier2019whole/species-trees-analysis"
WTREEQMC="$PROJDIR/software/TREE-QMC/tree-qmc"

cd $DATADIR

GTRE_FILE="UCE_minus_105_abayes_gene_trees_sorted_plus_bad105.tre"
PCSTREE="treeA-for-pcs"

# TREE-QMC w/ hybrid weighting scheme
MTHD="wtreeqmc_hybrid"
OUTF="${MTHD}_${PCSTREE}"
if [ ! -e $OUTF.tre ]; then
    $WTREEQMC --pcsonly $PCSTREE.tre -w h -n 0.333 -x 1 -d 0.333 --root galGal -i $GTRE_FILE -o $OUTF.tsv &> $OUTF.log
fi

# TREE-QMC w/ no weighting scheme
MTHD="wtreeqmc_none"
OUTF="${MTHD}_${PCSTREE}"
if [ ! -e $OUTF.tre ]; then
    $WTREEQMC --pcsonly $PCSTREE.tre -w n --root galGal -i $GTRE_FILE -o $OUTF.tsv &> $OUTF.log
fi

