#!/bin/bash

#exit


GTRES=( "introns_abayes_gene_trees_sorted.tre" \
        "nogal-introns_abayes_gene_trees_sorted.tre" \
	"nogal-UCEs_minus_bad105_sorted_plus_bad105_sorted_abayes_gene_trees.tre" \
	"UCEs_minus_bad105_sorted_plus_bad105_sorted_abayes_gene_trees.tre" \
	"nogal-UCEs_minus_bad105_abayes_gene_trees_sorted.tre" \
        "nogal-UCEs_minus_bad105_sorted_plus_bad105_sorted_without_galGal_and_tinGut_abayes_gene_trees.tre" \
        "UCEs_minus_bad105_abayes_gene_trees_sorted.tre" \
	"UCEs_minus_bad105_sorted_plus_bad105_sorted_without_galGal_and_tinGut_abayes_gene_trees.tre" )

for GTRE in ${GTRES[@]}; do
    echo "processing $GTRE"
    bash run_asteroid.sh $GTRE
    bash run_waster.sh $GTRE
    bash run_wastrid.sh $GTRE
    bash run_wtreeqmc.sh $GTRE
done

