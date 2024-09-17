#!/bin/bash

GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
DATADIR="$PROJDIR/data/cloutier2019whole/species-trees-analysis"
COMPARE="$PROJDIR/tools/compare_two_trees.py"

cd $DATADIR

STRE1="treeA.tre"
STRE2="treeA_w_hybrid_support_plus_bad105.tre"
python3 $COMPARE -t1 $STRE1 -t2 $STRE2

STRE1="treeA.tre"
STRE2="treeA_w_unweighted_support.tre"
python3 $COMPARE -t1 $STRE1 -t2 $STRE2

STRE1="treeA.tre"
STRE2="treeA_w_unweighted_support_plus_bad105.tre"
python3 $COMPARE -t1 $STRE1 -t2 $STRE2

STRE1="treeA.tre"
STRE2="treeA_w_unweighted_support_plus_bad105_without_galGal_tinGut.tre"
python3 $COMPARE -t1 $STRE1 -t2 $STRE2

STRE1="treeB.tre"
STRE2="treeB_w_hybrid_support_plus_bad105.tre"
python3 $COMPARE -t1 $STRE1 -t2 $STRE2

STRE1="treeB.tre"
STRE2="treeB_w_unweighted_support.tre"
python3 $COMPARE -t1 $STRE1 -t2 $STRE2

STRE1="treeB.tre"
STRE2="treeB_w_unweighted_support_plus_bad105.tre"
python3 $COMPARE -t1 $STRE1 -t2 $STRE2

STRE1="treeB.tre"
STRE2="treeB_w_unweighted_support_plus_bad105_without_galGal_tinGut.tre"
python3 $COMPARE -t1 $STRE1 -t2 $STRE2

