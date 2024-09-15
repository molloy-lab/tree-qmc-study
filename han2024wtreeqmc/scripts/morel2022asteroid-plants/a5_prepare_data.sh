#!/bin/bash

exit

# Define directories
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"

# Define software and utilities
RELABEL="$PROJDIR/tools/relabel_tree_by_map.py"

# Define input and output directories
INDIR="$GROUPDIR/group/data/morel2022asteroid/empirical/aa_plants_single"
OUTDIR="$PROJDIR/data/morel2022asteroid-plants"
GTRE="raxml-ng.LG+G.geneTree.newick"
STRE="speciesTree.newick"

if [ ! -e $OUTDIR/all-relabeled-$GTRE ]; then
    for GENE in `cat gene-family-ids.txt`; do
        cd $INDIR/families/$GENE
        pwd
        if [ ! -e relabeled-$GTRE ]; then
            python3 $RELABEL \
                -i gene_trees/$GTRE \
                -x mappings/treerecs_mapping.link \
                -o relabeled-$GTRE
        fi
    done

    cat $INDIR/families/family_*/relabeled-$GTRE > $OUTDIR/all-relabeled-$GTRE
fi

if [ ! -e $OUTDIR/reference-$STRE ]; then
    cp $INDIR/species_trees/$STRE $OUTDIR/reference-$STRE
    cp $INDIR/species_trees/concatenation-single.LG+G.speciesTree.newick $OUTDIR
fi

