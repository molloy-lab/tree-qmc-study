#!/bin/bash

exit

# ILS STUDY PARAMETERS
NTAXS=( 200 )                            # Number of taxa
HGHTS=( "10000000" "2000000" "500000" )  # Species tree height (number of generations)
RATES=( "0.0000001" "0.000001" )         # Speciation rate
REPLS=( $(seq -f "%02g" 1 50) )          # Replicates

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
 
echo "Submitting $MODL $REPL..."

sbatch \
    --nodelist="legacy[00-11,13-28,30]" \
    --job-name="d1.$MODL.$REPL" \
    --output="d1.$MODL.$REPL.%j.out" \
    --error="d1.$MODL.$REPL.%j.err" \
    --export=NTAX="$NTAX",HGHT="$HGHT",RATE="$RATE",REPL="$REPL" \
d1_drive.sbatch

#bash d1_compare_gene_trees.sh $NTAX $HGHT $RATE $REPL
            done
        done
    done
done

