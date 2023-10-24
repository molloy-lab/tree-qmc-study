#!/bin/bash

exit

MODLS=( $(cat model_list_sorted_focused.txt) )

for MODL in ${MODLS[@]}; do
    echo "Submitting $MODL..."
    sbatch \
        --job-name="h1.$MODL" \
        --output="h1.$MODL.%j.out" \
        --error="h1.$MODL.%j.err" \
        --export=MODL="$MODL" \
    h1_drive.sbatch
done
