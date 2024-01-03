#!/bin/bash

exit

# ILS STUDY PARAMETERS
NTAXS=( 200 )                            # Number of taxa
HGHTS=( "10000000" "2000000" "500000" )  # Species tree height (number of generations)
RATES=( "0.0000001" "0.000001" )         # Speciation rate

# SCALABILITY STUDY PARAMETERS
NTAXS=( 10 50 100 500 1000 )      # Number of taxa
HGHTS=( "2000000" )               # Species tree height (number of generations)
RATES=( "0.000001" )              # Speciation rate

# GENERAL PARAMETERS
REPLS=( $(seq -f "%02g" 1 50) )   # Replicates
NGENS=( 50 200 1000 )

for NTAX in ${NTAXS[@]}; do
    for HGHT in ${HGHTS[@]}; do
        for RATE in ${RATES[@]}; do
            MODL="model.$NTAX.$HGHT.$RATE"
            for REPL in ${REPLS[@]}; do
                for NGEN in ${NGENS[@]}; do
echo "Submitting $MODL $REPL $NGEN ..."
bash a1_copy_files.sh $NTAX $HGHT $RATE $REPL $NGEN
                done
            done
        done
    done
done

