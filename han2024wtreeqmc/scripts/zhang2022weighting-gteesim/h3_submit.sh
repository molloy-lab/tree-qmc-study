#!/bin/bash

exit

# STUDY PARAMETERS
REPLS=( $(seq -f "%02g" 1 50) )   # Replicates
NBPSS=( 200 400 800 1600 )
NGENS=( 50 200 500 1000 )

for NBPS in ${NBPSS[@]}; do
    for REPL in ${REPLS[@]}; do
        for NGEN in ${NGENS[@]}; do

MODL="model.$NBPS.$REPL.$NGEN"
echo "Submitting $MODL ..."

NAME="h3.$MODL"

sbatch \
    --job-name="$NAME" \
    --output="$NAME.%j.out" \
    --error="$NAME.%j.err" \
    --export=NBPS="$NBPS",REPL="$REPL",NGEN="$NGEN" \
h3_drive.sbatch
        done
    done
done

