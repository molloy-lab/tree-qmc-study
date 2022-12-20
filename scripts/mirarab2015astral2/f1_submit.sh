#!/bin/bash

exit

DATA="/fs/cbcb-lab/ekmolloy/ekmolloy/tree-qmc-study/data/mirarab2015astral2"


# ILS
NTAXS=( 200 )                             # Number of taxa
PSIZES=( "10000000" "2000000" "500000" )  # Population size
SRATES=( "0.0000001" "0.000001" )         # Speciation rate
REPLS=( $(seq  -f "%02g" 1 50) )          # Replicates
GTRE="estimatedgenetre"                   # Use estimated gene trees
#GTRE="truegenetrees"          # Use true gene trees
NGENS=( 250 1000 )                        # Number of gene trees

# SCALABILITY
#NTAXS=( 10 50 100 500 1000 )      # Number of taxa
#PSIZES=( "2000000" )              # Population size
#SRATES=( "0.000001" )             # Speciation rate
#REPLS=( $(seq  -f "%02g" 1 50) )  # Replicates
#GTRE="estimatedgenetre"           # Use estimated gene trees
##GTRE="truegenetrees"          # Use true gene trees
#NGENS=( 250 1000 )                # Number of gene trees

for NTAX in ${NTAXS[@]}; do
    for PSIZE in ${PSIZES[@]}; do
        for SRATE in ${SRATES[@]}; do
            MODL="model.$NTAX.$PSIZE.$SRATE"
                for NGEN in ${NGENS[@]}; do
                    for REPL in ${REPLS[@]}; do

echo "On $MODL $NGEN $REPL..."

FILEX="$DATA/$MODL/$REPL/treeqmc_n0_v1.0.0_estimatedgenetre_${NGEN}_node_info.csv"

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
        --job-name="f1.$MODL.$REPL.$GTRE.$NGEN" \
        --output="f1.$MODL.$REPL.$GTRE.$NGEN.%j.out" \
        --error="f1.$MODL.$REPL.$GTRE.$NGEN.%j.err" \
        --export=NTAX="$NTAX",PSIZE="$PSIZE",SRATE="$SRATE",REPL="$REPL",GTRE="$GTRE",NGEN="$NGEN" \
    f1_drive.sbatch
fi

                done
            done
        done
    done
done

