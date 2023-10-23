#!/bin/bash

exit

# Avian simulated data set

# (a) and (d)
#SCALS=( "0.5X" "1X" "2X" )
#NGENS=( 1000 )
#NBPSS=( "true" 500 )
#REPLS=( $(seq -f "R%g" 1 20) )

# (b) and (c)
#SCALS=( "1X" )
#NGENS=( 25 50 100 200 500 )
#NBPSS=( "true" 500 )
#REPLS=( $(seq -f "R%g" 1 20) )

for SCAL in ${SCALS[@]}; do
    for NGEN in ${NGENS[@]}; do
        for NBPS in ${NBPSS[@]}; do
            MODL="${SCAL}-${NGEN}-${NBPS}"
            for REPL in ${REPLS[@]}; do

            echo "Submitting $REPL..."
            sbatch \
                --job-name="g2.$MODL.$REPL" \
                --output="g2.$MODL.$REPL.%j.out" \
                --error="g2.$MODL.$REPL.%j.err" \
                --export=SCAL="$SCAL",NGEN="$NGEN",NBPS="$NBPS",REPL="$REPL" \
            g2_drive.sbatch

            done
        done
    done
done

