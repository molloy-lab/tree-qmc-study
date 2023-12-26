#!/bin/bash

exit

MODLS=( $(cat model_list_sorted_focused.txt) )

for MODL in ${MODLS[@]}; do
    echo "Submitting $MODL..."
    sbatch \
        --job-name="z1.$MODL" \
        --output="z1.$MODL.%j.out" \
        --error="z1.$MODL.%j.err" \
        --export=MODL="$MODL" \
    z1_drive.sbatch
done

