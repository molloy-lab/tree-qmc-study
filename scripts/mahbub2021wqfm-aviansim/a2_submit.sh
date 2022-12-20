#!/bin/bash

exit

# Avian simulated data set

SCALS=( "0.5X" "1X" "2X" )
REPLS=( $(seq -f "R%g" 1 20) )

for SCAL in ${SCALS[@]}; do
    for REPL in ${REPLS[@]}; do
        echo "Submitting $SCAL $REPL..."
        sbatch \
            --job-name="a2.$SCAL.$REPL" \
            --output="a2.$SCAL.$REPL.%j.out" \
            --error="a2.$SCAL.$REPL.%j.err" \
            --export=SCAL="$SCAL",REPL="$REPL" \
        a2_drive.sbatch
    done
done

