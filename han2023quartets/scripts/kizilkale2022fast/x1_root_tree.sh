#!/bin/bash

GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJECTDIR="$GROUPDIR/ekmolloy/species-to-tumor-study"

ROOTTREE="$PROJECTDIR/tools/root_tree_by_outgroup.py"
COMPARE="$PROJECTDIR/tools/compare_two_trees.py"

MYMODL=$1
DATDIR=$2

cd $DATDIR

pwd

MTHDS=( "treeqmcbip_v1.0.0_n2_wroot" \
        "fastme_v2.1.5_wroot" \
        "fastral_wroot" )

for MYMTHD in ${MTHDS[@]}; do
    FIN="$MYMTHD.cell_lineage_tree"
    FOUT="${MYMTHD}_wmuts"

    CTRE_ESTI_ROOTED="${MYMTHD}x.cell_lineage_tree.nwk"
    if [ -e $CTRE_ESTI ] && [ ! -e $CTRE_ESTI_ROOTED ]; then
        /opt/local/stow/Python3-3.8.1/bin/python3 $ROOTTREE \
            -i $CTRE_ESTI \
            -r "root" \
            -o  $CTRE_ESTI_ROOTED

        MYERRR=$(/opt/local/stow/Python3-3.8.1/bin/python3 $COMPARE -t1 $CTRE_ESTI -t2 $CTRE_ESTI_ROOTED)
        echo "$MYMODL,$MYMTHD,$MYERRR" > ${MYMTHD}x_cell_lineage_tree_error_debug.csv
    fi
done
