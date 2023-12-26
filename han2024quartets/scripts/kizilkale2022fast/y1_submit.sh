#!/bin/bash

exit

MODLS=( $(cat model_list_sorted_focused.txt) )

for MODL in ${MODLS[@]}; do
    echo "Submitting $MODL..."
    sbatch \
        --job-name="y1.$MODL" \
        --output="y1.$MODL.%j.out" \
        --error="y1.$MODL.%j.err" \
        --export=MODL="$MODL" \
    y1_drive.sbatch
done

