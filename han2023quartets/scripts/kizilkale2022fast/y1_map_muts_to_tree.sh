#!/bin/bash

GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJECTDIR="$GROUPDIR/ekmolloy/species-to-tumor-study"

MAPMUT2TRE="$PROJECTDIR/tools/map_mutations_to_tree.py"
BUILDTREE="$PROJECTDIR/tools/construct_tree_from_conflict_free_matrix.py"
COMPARE="$PROJECTDIR/tools/compare_two_trees.py"

MYMODL=$1
DATDIR=$2

cd $DATDIR

INPUT="noisy.CFMatrix"
CTRE_TRUE="ground.cell_lineage_tree.nwk"

pwd

MTHDS=( "treeqmcbip_v1.0.0_n2_wroot" \
        "scistree_v1.2.0.6" \
        "fastme_v2.1.5_wroot" \
        "fastral_wroot" )

for MYMTHD in ${MTHDS[@]}; do
    FIN="$MYMTHD.cell_lineage_tree"
    FOUT="${MYMTHD}_wmuts"

    if [ ! -e $FOUT.cell_lineage_tree.nwk ]; then    
        if [ -e $FIN.nwk ] && [ ! -e $FOUT.CFMatrix ]; then
            /opt/local/stow/Python3-3.8.1/bin/python3 $MAPMUT2TRE \
                -t $FIN.nwk \
                -m $INPUT \
                -o $FOUT.CFMatrix
        else
            echo "Either no file $FIN.nwk or else already file $FOUT.CFMatrix"
            exit
        fi

        if [ -e $FOUT.CFMatrix ]; then
            /opt/local/stow/Python3-3.8.1/bin/python3 $BUILDTREE \
                -i $FOUT.CFMatrix \
                -o $FOUT.cell_lineage_tree.nwk
        else
            echo "No file $FOUT.CFMatrix"
            exit
        fi

        #MYERRR=$(/opt/local/stow/Python3-3.8.1/bin/python3 $COMPARE -t1 $CTRE_TRUE -t2 $FIN.cell_lineage_tree.nwk)
        #echo "$MYMODL,${MYMTHD}_wmuts,$MYERRR" > ${MYMTHD}_cell_lineage_tree_error.csv

        MYERRR=$(/opt/local/stow/Python3-3.8.1/bin/python3 $COMPARE -t1 $CTRE_TRUE -t2 $FOUT.cell_lineage_tree.nwk)
        echo "$MYMODL,${MYMTHD}_wmuts,$MYERRR" > ${MYMTHD}_wmuts_cell_lineage_tree_error.csv
    fi
done

