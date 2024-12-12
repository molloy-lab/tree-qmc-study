#!/bin/bash

# Define directories and files
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
COMPARE="$PROJDIR/tools/compare_two_trees.py"
DATADIR="$PROJDIR/data/morel2022asteroid-plants"

cd $DATADIR


TREES1=( "apgweb_reference-speciesTree.newick" \
         "reference-speciesTree.newick" \
	 "concatenation-single.LG+G.speciesTree.newick" )

TREES2=( "concatenation-single.LG+G.speciesTree.newick" \
	 "treeqmc_wf_n2.tre" \
         "asteroid.bestTree.newick" \
         "wastrid_vanilla.tre" \
         "aster_v1.16.3.4.tre" )

for TREE1 in ${TREES1[@]}; do
    for TREE2 in ${TREES2[@]}; do
        echo "$TREE1 vs. $TREE2"
        python3 $COMPARE -t1 $TREE1 -t2 $TREE2
    done
done

exit

apgweb_reference-speciesTree.newick vs. concatenation-single.LG+G.speciesTree.newick
81,81,81,54,78,14,38
apgweb_reference-speciesTree.newick vs. treeqmc_wf_n2.tre
81,81,81,54,78,14,38
apgweb_reference-speciesTree.newick vs. asteroid.bestTree.newick
81,81,81,54,78,14,38
apgweb_reference-speciesTree.newick vs. wastrid_vanilla.tre
81,81,81,54,78,28,52
apgweb_reference-speciesTree.newick vs. aster_v1.16.3.4.tre
81,81,81,54,78,23,47

reference-speciesTree.newick vs. concatenation-single.LG+G.speciesTree.newick
81,81,81,56,78,7,29
reference-speciesTree.newick vs. treeqmc_wf_n2.tre
81,81,81,56,78,6,28
reference-speciesTree.newick vs. asteroid.bestTree.newick
81,81,81,56,78,6,28
reference-speciesTree.newick vs. wastrid_vanilla.tre
81,81,81,56,78,24,46
reference-speciesTree.newick vs. aster_v1.16.3.4.tre
81,81,81,56,78,17,39

concatenation-single.LG+G.speciesTree.newick vs. concatenation-single.LG+G.speciesTree.newick
81,81,81,78,78,0,0
concatenation-single.LG+G.speciesTree.newick vs. treeqmc_wf_n2.tre
81,81,81,78,78,15,15
concatenation-single.LG+G.speciesTree.newick vs. asteroid.bestTree.newick
81,81,81,78,78,14,14
concatenation-single.LG+G.speciesTree.newick vs. wastrid_vanilla.tre
81,81,81,78,78,29,29
concatenation-single.LG+G.speciesTree.newick vs. aster_v1.16.3.4.tre
81,81,81,78,78,20,20

