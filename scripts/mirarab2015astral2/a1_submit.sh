#!/bin/bash

exit

NTAXS=( 10 50 100 500 1000 )  # Number of taxa
NTAXS=( 1000 )
PSIZE="2000000"               # Population size
SRATE="0.000001"              # Speciation rate

for NTAX in ${NTAXS[@]}; do
    MODL="model.$NTAX.$PSIZE.$SRATE"
    sbatch \
        --job-name="a1.$MODL" \
        --output="a1.$MODL.%j.out" \
        --error="a1.$MODL.%j.err" \
        --export=NTAX="$NTAX",PSIZE="$PSIZE",SRATE="$SRATE" \
        a1_drive.sbatch
done

NTAX=200
PSIZES=( "10000000" "2000000" "500000" )
SRATES=( "0.0000001" "0.000001" )

for PSIZE in ${PSIZES[@]}; do
    for SRATE in ${SRATES[@]}; do
        MODL="model.$NTAX.$PSIZE.$SRATE"
        sbatch \
            --job-name="a1.$MODL" \
            --output="a1.$MODL.%j.out" \
            --error="a1.$MODL.%j.err" \
            --export=NTAX="$NTAX",PSIZE="$PSIZE",SRATE="$SRATE" \
            a1_drive.sbatch
    done
done

