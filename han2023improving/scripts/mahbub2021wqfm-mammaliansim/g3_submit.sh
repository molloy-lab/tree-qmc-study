#!/bin/bash

exit

# Mammalian simulated data set

# (a) varying ILS
#SCALS=( "noscale" "scale2d" "scale2u" )
#NGENS=( "200g" )
#NBPSS=( "500b" )

# (b) varying sequence length
#SCALS=("noscale" )
#NGENS=( "200g" )
#NBPSS=( "250b" "1000b" "1500b" "true" )

# (c) varying number of genes
#SCALS=( "noscale" )
#NGENS=( "25g" "50g" "100g" "400g" "800g" )
#NBPSS=( "500b" )

REPLS=( $(seq -f "R%g" 1 20) )

for SCAL in ${SCALS[@]}; do
    for NGEN in ${NGENS[@]}; do
        for NBPS in ${NBPSS[@]}; do
            MODL="${SCAL}.${NGEN}.${NBPS}"
            for REPL in ${REPLS[@]}; do

            echo "Submitting $REPL..."
            sbatch \
                --job-name="g3.$MODL.$REPL" \
                --output="g3.$MODL.$REPL.%j.out" \
                --error="g3.$MODL.$REPL.%j.err" \
                --export=SCAL="$SCAL",NGEN="$NGEN",NBPS="$NBPS",REPL="$REPL" \
            g3_drive.sbatch

            done
        done
    done
done

