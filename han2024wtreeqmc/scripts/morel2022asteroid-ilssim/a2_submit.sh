#!/bin/bash

#exit

MODLS=( $(cat model_list.txt) )

for MODL in ${MODLS[@]}; do
    echo "Submitting $MODL..."
    sbatch \
        --nodelist=$NODE \
        --job-name="a2.$MODL" \
        --output="a2.$MODL.%j.out" \
        --error="a2.$MODL.%j.err" \
        --export=MODL="$MODL" \
    a2_drive.sbatch
done

