#!/bin/bash

# Define software and utilities
GROUPDIR="/fs/cbcb-lab/ekmolloy"
HUNTRESS="$GROUPDIR/group/software/HUNTRESS-0.1.2.0/HUNTRESS.py"

PROJECTDIR="$GROUPDIR/ekmolloy/species-to-tumor-study"
BUILDTREE="$PROJECTDIR/tools/construct_tree_from_conflict_free_matrix.py"
COMPARE="$PROJECTDIR/tools/compare_two_trees.py"

MYMODL=$1
DATDIR=$2

MYMTHD="huntress_v0.1.2.0_default"
INPUT="noisy.CFMatrix"
CTRE_TRUE="ground.cell_lineage_tree.nwk"
CTRE_ESTI="$MYMTHD.cell_lineage_tree.nwk"
NTHREAD=16

cd $DATDIR
pwd

if [ ! -e $MYMTHD.CFMatrix ]; then
    MYTIME="$(time (/opt/local/stow/Python3-3.8.1/bin/python3 $HUNTRESS \
                    --i $INPUT \
                    --o $MYMTHD \
                    --t $NTHREAD \
                    &> $MYMTHD.log) 2>&1 1>/dev/null)"

    uname -a > ${MYMTHD}_node_info.csv
    MYNODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

    MYSECS="$(echo $MYTIME | awk '{print $2","$4","$6}')"
    echo "$MYMODL,$MYMTHD,$MYNODE,$NTHREAD,$MYSECS" > ${MYMTHD}_runtime.csv
fi

if [ -e $MYMTHD.CFMatrix ] && [ ! -e $CTRE_ESTI ]; then
        /opt/local/stow/Python3-3.8.1/bin/python3 $BUILDTREE \
            -i $MYMTHD.CFMatrix \
            -o $CTRE_ESTI

    MYERRR=$(/opt/local/stow/Python3-3.8.1/bin/python3 $COMPARE -t1 $CTRE_TRUE -t2 $CTRE_ESTI)
    echo "$MYMODL,$MYMTHD,$MYERRR" > ${MYMTHD}_cell_lineage_tree_error.csv
fi

