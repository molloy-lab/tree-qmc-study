#!/bin/bash

# Define software and utilities
GROUPDIR="/fs/cbcb-lab/ekmolloy"
SCISTREE="$GROUPDIR/group/software/ScisTree/ScisTree-ver1.2.0.6-src/scistree"

PROJECTDIR="$GROUPDIR/ekmolloy/species-to-tumor-study"
COMPARE="$PROJECTDIR/tools/compare_two_trees.py"

MYMODL=$1
DATDIR=$2

MYMTHD="scistree_v1.2.0.6"
INPUT="noisy.scistreemat"
CTRE_TRUE="ground.cell_lineage_tree.nwk"
CTRE_ESTI="$MYMTHD.cell_lineage_tree.nwk"
MTRE_ESTI="$MYMTHD.mutation_tree.nwk"
NTHREAD=1

cd $DATDIR
pwd

if [ ! -e $CTRE_ESTI ]; then
    MYTIME="$(time ($SCISTREE $INPUT -e -o $MYMTHD.mutation_tree.gml &> $MYMTHD.log) 2>&1 1>/dev/null)"

    uname -a > ${MYMTHD}_node_info.csv
    MYNODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

    MYSECS="$(echo $MYTIME | awk '{print $2","$4","$6}')"
    echo "$MYMODL,$MYMTHD,$MYNODE,$NTHREAD,$MYSECS" > ${MYMTHD}_runtime.csv

    TMP=$(grep "Mutation tree:" $MYMTHD.log | awk '{print $4}')
    echo "${TMP};" > $MTRE_ESTI

    TMP=$(grep "Constructed single cell phylogeny:" $MYMTHD.log | awk '{print $5}')
    echo "${TMP};" > $CTRE_ESTI

    MYERRR=$(/opt/local/stow/Python3-3.8.1/bin/python3 $COMPARE -t1 $CTRE_TRUE -t2 $CTRE_ESTI)
    echo "$MYMODL,$MYMTHD,$MYERRR" > ${MYMTHD}_cell_lineage_tree_error.csv
fi

