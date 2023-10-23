#!/bin/bash

exit

DATA="/fs/cbcb-lab/ekmolloy/ekmolloy/tree-qmc-study/data/mahbub2021wqfm/37-taxon-mammalian-simulated"

# Mammalian simulated data set

# (a) varying ILS
#SCALS=( "noscale" "scale2d" "scale2u" )
#NGENS=( "200g" )
#NBPSS=( "500b" )

# (b) varying sequence length
#SCALS=("noscale" )
#NGENS=( "200g" )
#NBPSS=( "250b" "1000b" "1500b" "true" )

# (c) varying number of genes
#SCALS=( "noscale" )
#NGENS=( "25g" "50g" "100g" "400g" "800g" )
#NBPSS=( "500b" )

REPLS=( $(seq -f "R%g" 1 20) )

for SCAL in ${SCALS[@]}; do
    for NGEN in ${NGENS[@]}; do
        for NBPS in ${NBPSS[@]}; do
            MODL="${SCAL}.${NGEN}.${NBPS}"
            for REPL in ${REPLS[@]}; do

echo "On $MODL $NGEN $REPL..."

FILEX="$DATA/$MODL/$REPL/treeqmc_n0_v1.0.0_node_info.csv"

if [ -e $FILEX ]; then
    NODE=$(cat $FILEX | sed 's/\./ /g' | awk '{print $2}')
    EXCLUDE="cbcbgpu00,cbcbgpu02,ibis[02-06,08-09,11-18,20-25],locust,redbud,spruce,tern[00-01],gum"
    if [ $NODE == "heron00" ]; then
        EXCLUDE="$EXCLUDE,heron[01-11]"
    elif [ $NODE == "heron01" ]; then
        EXCLUDE="$EXCLUDE,heron[00,02-11]"
    elif [ $NODE == "heron02" ]; then
        EXCLUDE="$EXCLUDE,heron[00-01,03-11]"
    elif [ $NODE == "heron03" ]; then
        EXCLUDE="$EXCLUDE,heron[00-02,04-11]"
    elif [ $NODE == "heron04" ]; then
        EXCLUDE="$EXCLUDE,heron[00-03,05-11]"
    elif [ $NODE == "heron05" ]; then
        EXCLUDE="$EXCLUDE,heron[00-04,06-11]"
    elif [ $NODE == "heron06" ]; then
        EXCLUDE="$EXCLUDE,heron[00-05,07-11]"
    elif [ $NODE == "heron07" ]; then
        EXCLUDE="$EXCLUDE,heron[00-06,08-11]"
    elif [ $NODE == "heron08" ]; then
        EXCLUDE="$EXCLUDE,heron[00-07,09-11]"
    elif [ $NODE == "heron09" ]; then
        EXCLUDE="$EXCLUDE,heron[00-08,10-11]"
    elif [ $NODE == "heron10" ]; then
        EXCLUDE="$EXCLUDE,heron[00-09,11]"
    elif [ $NODE == "heron11" ]; then
        EXCLUDE="$EXCLUDE,heron[00-10]"
    else
        echo "Unrecognized node info!"
        exit
    fi

    echo "... Found $NODE and running $EXCLUDE"

    sbatch \
        --exclude=$EXCLUDE \
        --job-name="c3.$MODL.$REPL" \
        --output="c3.$MODL.$REPL.%j.out" \
        --error="c3.$MODL.$REPL.%j.err" \
        --export=SCAL="$SCAL",NGEN="$NGEN",NBPS="$NBPS",REPL="$REPL" \
    c3_drive.sbatch
fi

            done
        done
    done
done

