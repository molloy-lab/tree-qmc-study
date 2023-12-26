#!/bin/bash

exit

# Define directories
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJECTDIR="$GROUPDIR/ekmolloy/species-to-tumor-study"

# Define software and utilities
BUILDTREE="$PROJECTDIR/tools/construct_tree_from_conflict_free_matrix.py"
ADDROOT="$PROJECTDIR/tools/add_root_to_matrix.py"
MATX2TREE="$PROJECTDIR/tools/convert_matrix_to_trees.py"
PREP4FASTME="$PROJECTDIR/tools/prepare_matrix_for_fastme.py"
PREP4SCISTREE="$PROJECTDIR/tools/prepare_matrix_for_scistree.py"

# Define data path
DATPATH="$PROJECTDIR/data/kizilkale2022fast"

for NSIM in `cat model_list_sorted_focused.txt`; do
    cd $DATPATH/$NSIM

    if [ -e ground.CFMatrix ]; then
        if [ ! -e ground.cell_lineage_tree.nwk ]; then
        /opt/local/stow/Python3-3.8.1/bin/python3 $BUILDTREE \
            -i ground.CFMatrix \
            -o ground.cell_lineage_tree.nwk
        fi

        if [ ! -e noisy_wroot.CFMatrix ]; then
        /opt/local/stow/Python3-3.8.1/bin/python3 $ADDROOT \
            -i noisy.CFMatrix \
            -o noisy_wroot.CFMatrix
        fi

        # Convert to newick strings
        if [ ! -e ground.nwk ]; then
        /opt/local/stow/Python3-3.8.1/bin/python3 $MATX2TREE \
            -i ground.CFMatrix \
            -o ground.nwk
        fi
        if [ ! -e noisy.nwk ]; then
        /opt/local/stow/Python3-3.8.1/bin/python3 $MATX2TREE \
            -i noisy.CFMatrix \
            -o noisy.nwk
        fi
        if [ ! -e noisy_wroot.nwk ]; then
        /opt/local/stow/Python3-3.8.1/bin/python3 $MATX2TREE \
            -i noisy_wroot.CFMatrix \
            -o noisy_wroot.nwk
        fi

        # Convert to FastME format
        if [ ! -e noisy.phy ]; then
        /opt/local/stow/Python3-3.8.1/bin/python3 $PREP4FASTME \
            -i noisy.CFMatrix \
            -o noisy.phy
        fi
        if [ ! -e noisy_wroot.phy ]; then
        /opt/local/stow/Python3-3.8.1/bin/python3 $PREP4FASTME \
            -i noisy_wroot.CFMatrix \
            -o noisy_wroot.phy
        fi

        # Convert to ScisTree format
        if [ ! -e noisy.scistreemat ]; then

        PARAMS=( $(echo $NSIM | sed 's/-/ /g' | sed 's/_/ /g') )
        END=$(echo "${#PARAMS[@]} - 1" | bc)
        for IND in `seq 0 2 $END`; do
            JND=$[IND+1]
            if [ ${PARAMS[$IND]} == "fp" ]; then
                FP=${PARAMS[$JND]}
            elif [ ${PARAMS[$IND]} == "fn" ]; then
                FN=${PARAMS[$JND]}
            fi
        done

        /opt/local/stow/Python3-3.8.1/bin/python3 $PREP4SCISTREE \
            -a $FP -b $FN \
            -i noisy.CFMatrix \
            -o noisy.scistreemat

        FP=""
        FN=""
        fi
    else
        echo "$DATPATH/$NSIM/ground.CFMatrix does not exist!"
    fi
done

