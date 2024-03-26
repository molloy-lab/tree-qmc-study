#!/bin/bash

exit

DATA="/fs/cbcb-lab/ekmolloy/ekmolloy/tree-qmc-study/han2024wtreeqmc/data/mirarab2015astral2-extsim/"

# ILS STUDY PARAMETERS
NTAXS=( 200 )                            # Number of taxa
HGHTS=( "10000000" "2000000" "500000" )  # Species tree height (number of generations)
RATES=( "0.0000001" "0.000001" )         # Speciation rate

# SCALABILITY STUDY PARAMETERS
NTAXS=( 10 50 100 500 1000 )      # Number of taxa
HGHTS=( "2000000" )               # Species tree height (number of generations)
RATES=( "0.000001" )              # Speciation rate

# GENERAL PARAMETERS
REPLS=( $(seq -f "%02g" 1 50) )   # Replicates
SUPPS=( "abayes" )
NGENS=( 1000 200 50 )

for NTAX in ${NTAXS[@]}; do
    for HGHT in ${HGHTS[@]}; do
        for RATE in ${RATES[@]}; do
            for REPL in ${REPLS[@]}; do
                for SUPP in ${SUPPS[@]}; do
    		    for NGEN in ${NGENS[@]}; do

MODL="model.$NTAX.$HGHT.$RATE.$REPL.$SUPP.$NGEN"
echo "Submitting $MODL ..."

FILEX="$DATA/model.$NTAX.$HGHT.$RATE/$REPL/wastrid_s_${SUPP}_${NGEN}gen_node_info.csv"

if [ -e $FILEX ]; then
    NODE=$(cat $FILEX | sed 's/\./ /g' | awk '{print $2}')
    echo "  using node $NODE"
    sbatch \
        --nodelist=$NODE \
        --job-name="etoh1.$MODL" \
        --output="etoh1.$MODL.%j.out" \
        --error="etoh1.$MODL.%j.err" \
        --export=NTAX="$NTAX",HGHT="$HGHT",RATE="$RATE",REPL="$REPL",SUPP="$SUPP",NGEN="$NGEN" \
    etoh1_drive.sbatch
else
    sbatch \
        --job-name="etoh1.$MODL" \
        --output="etoh1.$MODL.%j.out" \
        --error="etoh1.$MODL.%j.err" \
        --export=NTAX="$NTAX",HGHT="$HGHT",RATE="$RATE",REPL="$REPL",SUPP="$SUPP",NGEN="$NGEN" \
    etoh1_drive.sbatch
fi
                    done
                done
            done
        done
    done
done

