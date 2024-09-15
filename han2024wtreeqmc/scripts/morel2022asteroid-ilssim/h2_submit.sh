#!/bin/bash

exit

MODLS=( $(cat model_list.txt) )
MODLS=( $(cat model_list_s25.txt) )
MODLS=( $(cat model_list_s75.txt) )
MODLS=( $(cat model_list_s100.txt) )
MODLS=( $(cat model_list_s125.txt) )
MODLS=( $(cat model_list_s150.txt) )
MODLS=( $(head -n210 model_list_s50.txt) )
MODLS=( $(head -n420 model_list_s50.txt | tail -n210) )
MODLS=( $(head -n630 model_list_s50.txt | tail -n210) )
MODLS=( $(head -n840 model_list_s50.txt | tail -n210) )
MODLS=( $(head -n1050 model_list_s50.txt | tail -n210) )

for MODL in ${MODLS[@]}; do
    echo "Submitting $MODL..."

    sbatch \
        --nodelist=$NODE \
        --job-name="h2.$MODL" \
        --output="h2.$MODL.%j.out" \
        --error="h2.$MODL.%j.err" \
        --export=MODL="$MODL" \
    h2_drive.sbatch

done

