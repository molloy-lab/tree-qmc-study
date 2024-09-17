#!/bin/bash

TYPE=$1  # Data type
GENE=$2  # Gene ID

# Define input files
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
DATADIR="$PROJDIR/data/cloutier2019whole"

# Do work
cd $DATADIR/RAxML_bestML_gene_trees/$TYPE

#pwd

# Check inputs
if [ ! -e $GENE.tre.abayes ]; then
    echo "$GENE.tre.abayes does not exist!"
    exit 1
fi

# Summarize branch support after mapping [0.333, 1] to [0, 1]
if [ ! -e $GENE.csv ] || [ -z $(cat $GENE.csv) ]; then
	echo $GENE
	python3 $PROJDIR/tools/get_branch_support_summary.py \
		-i $GENE.tre.abayes \
		-p "$TYPE,$GENE" \
		> $GENE.csv
fi

