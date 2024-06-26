#!/bin/bash

exit

# SCALABILITY STUDY PARAMETERS
NTAXS=( 10 50 100 500 1000 )      # Number of taxa
HGHTS=( "2000000" )               # Species tree height (number of generations)
RATES=( "0.000001" )              # Speciation rate
REPLS=( $(seq -f "%02g" 1 50) )   # Replicates

for NTAX in ${NTAXS[@]}; do
    for HGHT in ${HGHTS[@]}; do
        for RATE in ${RATES[@]}; do
            MODL="model.$NTAX.$HGHT.$RATE"
            for REPL in ${REPLS[@]}; do

echo "Running $MODL $REPL..."
./c1_organize_gene_trees.sh $MODL $REPL

            done
        done
    done
done

