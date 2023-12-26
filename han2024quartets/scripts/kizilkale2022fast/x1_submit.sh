#!/bin/bash

exit

MODLS=( $(cat model_list_sorted_focused.txt) )

for MODL in ${MODLS[@]}; do
    echo "Submitting $MODL..."
    sbatch \
        --job-name="x1.$MODL" \
        --output="x1.$MODL.%j.out" \
        --error="x1.$MODL.%j.err" \
        --export=MODL="$MODL" \
    x1_drive.sbatch
done

