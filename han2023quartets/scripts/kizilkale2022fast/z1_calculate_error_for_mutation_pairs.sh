#!/bin/bash

GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJECTDIR="$GROUPDIR/ekmolloy/species-to-tumor-study"

CALCERROR="$PROJECTDIR/tools/calculate_mutation_error.py "

MYMODL=$1
DATDIR=$2

cd $DATDIR

TRUE="ground.CFMatrix"

pwd

MTHDS=( "huntress_v0.1.2.0_default" \
        "treeqmcbip_v1.0.0_n2_wrootx_wmuts" \
        "scistree_v1.2.0.6_wmuts" \
        "fastme_v2.1.5_wrootx_wmuts" \
        "fastral_wrootx_wmuts" )

MTHDS=( "aster_v1.10.2.1_wrootx_wmuts" \
        "treeqmcbip_v1.0.0_n0_wrootx_wmuts" \
        "treeqmcbip_v1.0.0_n1_wrootx_wmuts" )

for MYMTHD in ${MTHDS[@]}; do
    ESTI="$MYMTHD.CFMatrix"
    FOUT="${MYMTHD}_mut_pair_error.csv"
    if [ -e $ESTI ] && [ ! -e $FOUT ]; then    
        MYERRR=$(/opt/local/stow/Python3-3.8.1/bin/python3 $CALCERROR -t $TRUE -e $ESTI)
        echo "$MYMODL,$MYMTHD,$MYERRR" > $FOUT
    fi
done

