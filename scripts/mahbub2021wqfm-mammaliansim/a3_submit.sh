#!/bin/bash

exit

# Mammalian simulated data set

NBPSS=( "250b" "500b" "1000b" "1500b" )
REPLS=( $(seq -f "R%g" 1 20) )

for NBPS in ${NBPSS[@]}; do
    for REPL in ${REPLS[@]}; do
        echo "Submitting $NBPS $REPL..."
        sbatch \
            --job-name="a3.$NBPS.$REPL" \
            --output="a3.$NBPS.$REPL.%j.out" \
            --error="a3.$NBPS.$REPL.%j.err" \
            --export=NBPS="$NBPS",REPL="$REPL" \
        a3_drive.sbatch
    done
done

