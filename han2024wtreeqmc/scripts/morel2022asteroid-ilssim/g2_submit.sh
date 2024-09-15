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

MODLS=( "ssim_veryhighmiss_s50_f1000_sites100_GTR_bl1.0_d0.0_l0.0_t0.0_gc0.0_p0.0_pop50000000_ms0.75_mf0.75_seed3047" \
        "ssim_veryhighmiss_s50_f1000_sites100_GTR_bl1.0_d0.0_l0.0_t0.0_gc0.0_p0.0_pop50000000_ms0.75_mf0.75_seed3018" \
	"ssim_veryhighmiss_s50_f1000_sites100_GTR_bl1.0_d0.0_l0.0_t0.0_gc0.0_p0.0_pop10_ms0.6_mf0.6_seed3037" )

for MODL in ${MODLS[@]}; do
    echo "Submitting $MODL..."

    sbatch \
        --nodelist=$NODE \
        --job-name="g2.$MODL" \
        --output="g2.$MODL.%j.out" \
        --error="g2.$MODL.%j.err" \
        --export=MODL="$MODL" \
    g2_drive.sbatch

done

