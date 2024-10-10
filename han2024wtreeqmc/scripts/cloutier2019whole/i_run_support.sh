#!/bin/bash

# Define directories
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
GTREDIR="$PROJDIR/data/cloutier2019whole/gene-trees"
DATADIR="$PROJDIR/data/cloutier2019whole/species-trees-analysis"
ASTERH="$GROUPDIR/group/software-compiled-on-EPYC-7313/ASTER/bin/astral-hybrid"
#--mode	Integer	1	1: hybrid weighting, 2: support only, 3: length only, 4: unweighted
#--scoring	Preset		Scoring the full species tree file after `-c` without exploring other topologies (`-r 1 -s 0`)

cd $DATADIR

for STRE in "treeA" "treeB"; do

# Tree w/ hybrid weighting scheme (full input)
GTRE="UCEs_minus_bad105_sorted_plus_bad105_sorted_abayes_gene_trees"
$ASTERH --mode 1 --bayes \
        -t 16 -u 2 -r 1 -s 0 \
        --scoring --root galGal \
        -c $STRE.tre \
        -i $GTREDIR/${GTRE}.tre \
        -o ${STRE}_w_hybrid_support_plus_bad105.tre \
        &> ${STRE}_w_hybrid_support_plus_bad105.log

# Tree w/ unweighted scheme (input after filtering taxa with homology errors)
GTRE="UCEs_minus_bad105_sorted_plus_bad105_sorted_without_galGal_and_tinGut_abayes_gene_trees"
$ASTERH --mode 4 \
        -t 16 -u 2 -r 1 -s 0 \
        --scoring --root galGal \
        -c $STRE.tre \
        -i $GTREDIR/${GTRE}.tre \
        -o ${STRE}_w_unweighted_support_plus_bad105_without_galGal_tinGut.tre \
        &> ${STRE}_w_unweighted_support_plus_bad105_without_galGal_tinGut.log

# Tree w/ unweighted scheme (input after filtering genes with homology errors)
GTRE="UCEs_minus_bad105_abayes_gene_trees_sorted"
$ASTERH --mode 4 \
        -t 16 -u 2 -r 1 -s 0 \
        --scoring --root galGal \
        -c $STRE.tre \
        -i $GTREDIR/${GTRE}.tre \
        -o ${STRE}_w_unweighted_support.tre \
        &> ${STRE}_w_unweighted_support.log

# Tree w/ unweighted weighting scheme (full input)
GTRE="UCEs_minus_bad105_sorted_plus_bad105_sorted_abayes_gene_trees"
$ASTERH --mode 4 \
        -t 16 -u 2 -r 1 -s 0 \
        --scoring --root galGal \
        -c $STRE.tre \
        -i $GTREDIR/${GTRE}.tre \
        -o ${STRE}_w_unweighted_support_plus_bad105.tre \
        &> ${STRE}_w_unweighted_support_plus_bad105.log

done

#treeA_w_hybrid_support_plus_bad105.log:Score: 2326403.019
#treeA_w_unweighted_support.log:Score: 3031323
#treeA_w_unweighted_support_plus_bad105.log:Score: 3117457
#treeA_w_unweighted_support_plus_bad105_without_galGal_tinGut.log:Score: 3083056
#treeB_w_hybrid_support_plus_bad105.log:Score: 2325430.563
#treeB_w_unweighted_support.log:Score: 3029054
#treeB_w_unweighted_support_plus_bad105.log:Score: 3119240
#treeB_w_unweighted_support_plus_bad105_without_galGal_tinGut.log:Score: 3080549

