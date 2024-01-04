#!/bin/bash

#exit

# STUDY PARAMETERS
REPLS=( $(seq -f "%02g" 1 50) )   # Replicates
NBPSS=( 200 400 800 1600 )
SUPPS=( "bs" "abayes" )
NGENS=( 50 200 500 1000 )

for REPL in ${REPLS[@]}; do
    for SUPP in ${SUPPS[@]}; do
        for NBPS in ${NBPSS[@]}; do
            for NGEN in ${NGENS[@]}; do

MODL="$REPL.$NBPS.$SUPP.$NGEN"
echo "Submitting $MODL ..."

sbatch \
    --job-name="btog3.$MODL" \
    --output="btog3.$MODL.%j.out" \
    --error="btog3.$MODL.%j.err" \
    --export=REPL="$REPL",NBPS="$NBPS",SUPP="$SUPP",NGEN="$NGEN" \
btog3_drive.sbatch

#bash btog3_test.sh $REPL $NBPS $SUPP $NGEN
            done
        done
    done
done

