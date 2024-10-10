#!/bin/bash

exit

GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
DATADIR="$PROJDIR/data/cloutier2019whole/RAxML_bestML_gene_trees"
OUTDIR="$PROJDIR/data/cloutier2019whole/gene-trees"
RESTRICT="$PROJDIR/tools/restrict_taxa.py"

cd $DATADIR

# Handle CNEEs
INPUT="CNEEs_abayes_gene_trees_sorted.tre"
if [ ! -e $INPUT ]; then
    mv CNEEs/$INPUT $OUTDIR
fi

# Handle introns
INPUT="introns_abayes_gene_trees_sorted.tre"
if [ ! -e $INPUT ]; then
    mv introns/$INPUT $OUTDIR
fi

# Handle good UCEs
INPUT="UCEs_minus_bad105_abayes_gene_trees_sorted.tre"
if [ ! -e $INPUT ]; then
    mv UCEs/$INPUT $OUTDIR
fi

# Remove impacted taxa from bad UCEs
cd UCEs/exclude
INPUT="bad105_abayes_gene_trees_sorted"
if [ ! -e ${INPUT}_without_galGal_tinGut.tre ]; then
    python3 $RESTRICT -i $INPUT.tre -r "galGal,tinGut" > ${INPUT}_without_galGal_tinGut.tre
fi
cd ../../

# Combine good and bad UCEs
OUTPUT="UCEs_minus_bad105_sorted_plus_bad105_sorted_abayes_gene_trees.tre"
if [ ! -e $OUTPUT ]; then
    cat UCEs_minus_bad105_abayes_gene_trees_sorted.tre UCEs/exclude/bad105_abayes_gene_trees_sorted.tre > $OUTDIR/$OUTPUT
fi

OUTPUT="UCEs_minus_bad105_sorted_plus_bad105_sorted_without_galGal_and_tinGut_abayes_gene_trees.tre"
if [ ! -e $OUTPUT ]; then
    cat UCEs_minus_bad105_abayes_gene_trees_sorted.tre UCEs/exclude/bad105_abayes_gene_trees_sorted_without_galGal_tinGut.tre > $OUTDIR/$OUTPUT
fi

