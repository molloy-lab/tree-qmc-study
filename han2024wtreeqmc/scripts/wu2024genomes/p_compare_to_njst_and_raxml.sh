#!/bin/bash

GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="/fs/cbcb-lab/ekmolloy/ekmolloy/tree-qmc-study/han2024wtreeqmc"
OUTDIR="$PROJDIR/data/wu2024genomes/species-trees"
COMPARE="$PROJDIR/tools/compare_two_trees.py"

cd $OUTDIR

NJST="wu2024_rooted_njst_threeMarkers_70%_outliers_contree.tre"
RAXML="wu2024_rooted_raxml_threeMarkers_70%_outliers_contree.tre"

TREEA="rooted_treeA_wastrid_support_none_nonefilter.tre"
TREEB="rooted_treeB_wtreeqmc_hybrid_length.tre"
TREEC="rooted_treeD_waster_hybrid.tre"
TREED="rooted_treeC_wtreeqmc_support_none_nonefilter.tre"
TREEE="rooted_treeG_waster_nonefilter.tre"

echo "b) ASTRID and Asteroid (after treeshrink)"
python3 $COMPARE -t1 $NJST -t2 $TREEA
#125,125,125,122,122,8,8
python3 $COMPARE -t1 $RAXML -t2 $TREEA
#125,125,125,122,122,16,16

echo "c) weighted TREEQMC (hybrid)"
python3 $COMPARE -t1 $NJST -t2 $TREEB
# 125,125,125,122,122,8,8
python3 $COMPARE -t1 $RAXML -t2 $TREEB
# 125,125,125,122,122,13,13

echo "d) weighted ASTER (hybrid)"
python3 $COMPARE -t1 $NJST -t2 $TREEC
# 125,125,125,122,122,16,16
python3 $COMPARE -t1 $RAXML -t2 $TREEC
# 125,125,125,122,122,18,18

echo "e) TREE-QMC"
python3 $COMPARE -t1 $NJST -t2 $TREED
# 125,125,125,122,122,11,11
python3 $COMPARE -t1 $RAXML -t2 $TREED
# 125,125,125,122,122,13,13

echo "e) ASTER (after TreeShrink)"
python3 $COMPARE -t1 $NJST -t2 $TREEE
# 125,125,125,122,122,14,14
python3 $COMPARE -t1 $RAXML -t2 $TREEE
# 125,125,125,122,122,16,16

