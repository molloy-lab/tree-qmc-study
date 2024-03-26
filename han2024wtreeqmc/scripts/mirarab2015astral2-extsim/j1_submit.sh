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
REPLS=( $(seq -f "%02g" 1 50) )  # Replicates
SUPPS=( "abayes" )
NGENS=( 1000 )

for NTAX in ${NTAXS[@]}; do
    for HGHT in ${HGHTS[@]}; do
        for RATE in ${RATES[@]}; do
            for REPL in ${REPLS[@]}; do
                for SUPP in ${SUPPS[@]}; do
                    for NGEN in ${NGENS[@]}; do

MODL="model.$NTAX.$HGHT.$RATE.$REPL.$SUPP.$NGEN"
echo "Submitting $MODL ..."

NAME="j1.$MODL"

sbatch \
    --job-name="$NAME" \
    --output="$NAME.%j.out" \
    --error="$NAME.%j.err" \
    --export=NTAX="$NTAX",HGHT="$HGHT",RATE="$RATE",REPL="$REPL",SUPP="$SUPP",NGEN="$NGEN" \
j1_drive.sbatch

#exit
                    done
                done
            done
        done
    done
done

