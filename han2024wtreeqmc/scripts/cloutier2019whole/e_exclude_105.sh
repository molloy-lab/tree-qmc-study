#!/bin/bash

exit

GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
DATADIR="$PROJDIR/data/cloutier2019whole"
EXCLUDE="$PROJDIR/scripts/cloutier2019whole/list_of_105_UCEs_to_exclude.txt"

cd $DATADIR/RAxML_bestML_gene_trees/UCEs
mkdir exclude
for GENE in `cat $EXCLUDE`; do
    mv $GENE.tre.abayes exclude.$GENE.tre.abayes
    mv $GENE.csv exclude.$GENE.csv
    mv exclude.* exclude
done

cd $DATADIR/alignments/UCEs
mkdir exclude
for GENE in `cat $EXCLUDE`; do
    mv $GENE.fasta exclude.$GENE.fasta
    mv $GENE.csv exclude.$GENE.csv
    mv exclude.* exclude
done

