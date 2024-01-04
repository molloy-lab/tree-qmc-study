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
REPLS=( $(seq -f "%02g" 1 25) )   # Replicates
SUPPS=( "sh" "abayes" )
NGENS=( 1000 200 50 )

for NTAX in ${NTAXS[@]}; do
    for HGHT in ${HGHTS[@]}; do
        for RATE in ${RATES[@]}; do
            for REPL in ${REPLS[@]}; do
                for SUPP in ${SUPPS[@]}; do
                    for NGEN in ${NGENS[@]}; do

MODL="model.$NTAX.$HGHT.$RATE.$REPL.$SUPP.$NGEN"
echo "Submitting $MODL ..."

sbatch \
    --job-name="h1.$MODL" \
    --output="h1.$MODL.%j.out" \
    --error="h1.$MODL.%j.err" \
    --export=NTAX="$NTAX",HGHT="$HGHT",RATE="$RATE",REPL="$REPL",SUPP="$SUPP",NGEN="$NGEN" \
h1_drive.sbatch

                    done
                done
            done
        done
    done
done

