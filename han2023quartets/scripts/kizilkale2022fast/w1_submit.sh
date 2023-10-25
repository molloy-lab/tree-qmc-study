#!/bin/bash

exit

MODLS=( $(cat model_list_sorted_focused.txt) )

for MODL in ${MODLS[@]}; do
    echo "Submitting $MODL..."
    sbatch \
        --job-name="w1.$MODL" \
        --output="w1.$MODL.%j.out" \
        --error="w1.$MODL.%j.err" \
        --export=MODL="$MODL" \
    w1_drive.sbatch
done

