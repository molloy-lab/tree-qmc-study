#!/bin/bash

exit

# STUDY PARAMETERS
REPLS=( $(seq -f "%02g" 1 50) )   # Replicates
NBPSS=( 200 400 800 1600 )

for REPL in ${REPLS[@]}; do
    for NBPS in ${NBPSS[@]}; do

MODL="$REPL.$NBPS"
echo "Submitting $MODL ..."

    sbatch \
        --nodelist="legacy[00-11,13-28,30]" \
        --job-name="a3.$MODL" \
        --output="a3.$MODL.%j.out" \
        --error="a3.$MODL.%j.err" \
        --export=REPL="$REPL",NBPS="$NBPS" \
    a3_drive.sbatch

    done
done

