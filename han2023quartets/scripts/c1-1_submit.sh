#!/bin/bash

exit

DATA="/fs/cbcb-lab/ekmolloy/ekmolloy/species-to-tumor-study/data/kizilkale2022fast"
MODLS=( $(cat model_list_sorted_focused.txt) )
NODES=( $(seq -f "heron%02g" 3 11) )
NNODES=${#NODES[@]}

INDEX=0
for MODL in ${MODLS[@]}; do
    echo "Submitting $MODL..."

    if [ $INDEX -eq $NNODES ]; then
        INDEX=0
    fi

    FILEX="$DATA/$MODL/huntress_v0.1.2.0_default_node_info.csv"
    if [ ! -e $FILEX ]; then
        NODE="${NODES[$INDEX]}"
        echo "Should be submitted to $NODE"

        sbatch \
            --nodelist=$NODE \
            --job-name="c1-1.$MODL" \
            --output="c1-1.$MODL.%j.out" \
            --error="c1-1.$MODL.%j.err" \
            --export=MODL="$MODL" \
        c1-1_drive.sbatch

        INDEX=$[INDEX+1]
    else
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

    #sbatch \
    #    --exclude=$EXCLUDE \
    #    --job-name="c1-1.$MODL" \
    #    --output="c1-1.$MODL.%j.out" \
    #    --error="c1-1.$MODL.%j.err" \
    #    --export=MODL="$MODL" \
    #c1-1_drive.sbatch

    sbatch \
        --nodelist=$NODE \
        --job-name="c1-1.$MODL" \
        --output="c1-1.$MODL.%j.out" \
        --error="c1-1.$MODL.%j.err" \
        --export=MODL="$MODL" \
    c1-1_drive.sbatch

    fi
done
