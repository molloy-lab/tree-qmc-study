#!/bin/bash

GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
DATADIR="$PROJDIR/data/cloutier2019whole/species-trees-analysis"
COMPARE="$PROJDIR/tools/compare_two_trees.py"

cd $DATADIR

# Check - TREE-QMC trees
T1="wtreeqmc_hybrid_n2_UCE_minus_105_abayes_gene_trees_sorted_plus_bad105.tre"
T2="wtreeqmc_length_n2_UCE_minus_105_abayes_gene_trees_sorted_plus_bad105.tre"
T3="wtreeqmc_support_n2_UCE_minus_105_abayes_gene_trees_sorted_plus_bad105.tre"
T4="wtreeqmc_none_n2_UCE_minus_105_abayes_gene_trees_sorted_plus_bad105.tre"
T5="wtreeqmc_none_n2_UCE_minus_105_abayes_gene_trees_sorted_plus_bad105_without_galGal_tinGut.tre"
T6="wtreeqmc_none_n2_UCE_minus_105_abayes_gene_trees_sorted.tre"

#python3 $COMPARE -t1 $T1 -t2 $T2
# 15,15,15,12,12,0,0

#python3 $COMPARE -t1 $T1 -t2 $T3
#15,15,15,12,12,1,1

#python3 $COMPARE -t1 $T1 -t2 $T4
#15,15,15,12,12,1,1

#python3 $COMPARE -t1 $T1 -t2 $T5
#15,15,15,12,12,0,0

#python3 $COMPARE -t1 $T1 -t2 $T6
#15,15,15,12,12,0,0

#python3 $COMPARE -t1 $T3 -t2 $T4
#15,15,15,12,12,0,0

TA=$T1 # TREE-QMC (hybrid weighting, length weighting, taxa filtering, gene filtering)
       # ASTRAL   (hybrid weighting, length weighting, taxa filtering, gene filtering)
TB=$T3 # TREE-QMC (no weighting, support weighting); 
       # ASTRAL   (no weighting, support weighting); 
       # ASTRID   (no weighting, support weighting, length weighting, taxa filtering, gene filtering)

# Check ASTER trees
T7="waster_hybrid_UCE_minus_105_abayes_gene_trees_sorted_plus_bad105.tre"
T8="waster_length_UCE_minus_105_abayes_gene_trees_sorted_plus_bad105.tre"
T9="waster_support_UCE_minus_105_abayes_gene_trees_sorted_plus_bad105.tre"
T10="waster_none_UCE_minus_105_abayes_gene_trees_sorted_plus_bad105.tre"
T11="waster_none_UCE_minus_105_abayes_gene_trees_sorted_plus_bad105_without_galGal_tinGut.tre"
T12="waster_none_UCE_minus_105_abayes_gene_trees_sorted.tre"

#python3 $COMPARE -t1 $TA -t2 $T7
#15,15,15,12,12,0,0

#python3 $COMPARE -t1 $TA -t2 $T8
#15,15,15,12,12,0,0

#python3 $COMPARE -t1 $TB -t2 $T9
#15,15,15,12,12,0,0

#python3 $COMPARE -t1 $TB -t2 $T10
#15,15,15,12,12,0,0

#python3 $COMPARE -t1 $TA -t2 $T11
#15,15,15,12,12,0,0

#python3 $COMPARE -t1 $TA -t2 $T12
#15,15,15,12,12,0,0

T13="wastrid_length_UCE_minus_105_abayes_gene_trees_sorted_plus_bad105.tre"
T14="wastrid_support_UCE_minus_105_abayes_gene_trees_sorted_plus_bad105.tre"
T15="wastrid_none_UCE_minus_105_abayes_gene_trees_sorted_plus_bad105.tre"
T16="wastrid_none_UCE_minus_105_abayes_gene_trees_sorted_plus_bad105_without_galGal_tinGut.tre"
T17="wastrid_none_UCE_minus_105_abayes_gene_trees_sorted.tre"

#python3 $COMPARE -t1 $TB -t2 $T13
#15,15,15,12,12,0,0

#python3 $COMPARE -t1 $TB -t2 $T14
#15,15,15,12,12,0,0

#python3 $COMPARE -t1 $TB -t2 $T15
#15,15,15,12,12,0,0

#python3 $COMPARE -t1 $TB -t2 $T16
#15,15,15,12,12,0,0

#python3 $COMPARE -t1 $TB -t2 $T16
#15,15,15,12,12,0,0

T18="asteroid_UCE_minus_105_abayes_gene_trees_sorted_plus_bad105.tre.bestTree.newick"
T19="asteroid_UCE_minus_105_abayes_gene_trees_sorted_plus_bad105_without_galGal_tinGut.tre.bestTree.newick"
T20="asteroid_UCE_minus_105_abayes_gene_trees_sorted.tre.bestTree.newick"

#python3 $COMPARE -t1 $TB -t2 $T18
#15,15,15,12,12,0,0

#python3 $COMPARE -t1 $TB -t2 $T19
#15,15,15,12,12,0,0

#python3 $COMPARE -t1 $TB -t2 $T20
#15,15,15,12,12,0,0

# Lastly, compare to published RAxML and ExaML trees
#python3 $COMPARE -t1 cloutier2019_UCE_examl.tre -t2 simmons2022_UCE_minus_105_raxml.tre 
#15,15,15,12,12,0,0

#python3 $COMPARE -t1 $TA -t2 cloutier2019_UCE_examl.tre
#15,15,15,12,12,1,1

#python3 $COMPARE -t1 $TB -t2 cloutier2019_UCE_examl.tre
#15,15,15,12,12,2,2

