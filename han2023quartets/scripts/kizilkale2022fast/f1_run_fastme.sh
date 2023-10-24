#!/bin/bash

# Define software and utilities
GROUPDIR="/fs/cbcb-lab/ekmolloy"
FASTME="$GROUPDIR/group/software/fastme-2.1.5/binaries/fastme-2.1.5-linux64-omp"

PROJECTDIR="$GROUPDIR/ekmolloy/species-to-tumor-study"
ROOTTREE="$PROJECTDIR/tools/root_tree.py"
COMPARE="$PROJECTDIR/tools/compare_two_trees.py"

MYMODL=$1
DATDIR=$2

#MYMTHD="fastme_v2.1.5"
MYMTHD="fastme_v2.1.5_wroot"

#INPUT="noisy.phy"
INPUT="noisy_wroot.phy"
CTRE_TRUE="ground.cell_lineage_tree.nwk"
CTRE_ESTI="$MYMTHD.cell_lineage_tree.nwk"
NTHREAD=16

cd $DATDIR
pwd

if [ ! -e $CTRE_ESTI ]; then
    MYTIME="$(time ($FASTME -i $INPUT -dp -mB --nni=B --spr -T $NTHREAD -o $CTRE_ESTI &> $MYMTHD.log) 2>&1 1>/dev/null)"

    uname -a > ${MYMTHD}_node_info.csv
    MYNODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

    MYSECS="$(echo $MYTIME | awk '{print $2","$4","$6}')"
    echo "$MYMODL,$MYMTHD,$MYNODE,$NTHREAD,$MYSECS" > ${MYMTHD}_runtime.csv

    MYERRR=$(/opt/local/stow/Python3-3.8.1/bin/python3 $COMPARE -t1 $CTRE_TRUE -t2 $CTRE_ESTI)
    echo "$MYMODL,$MYMTHD,$MYERRR" > ${MYMTHD}_cell_lineage_tree_error.csv
fi

CTRE_ESTI_ROOTED="${MYMTHD}x.cell_lineage_tree.nwk"
if [ -e $CTRE_ESTI ] && [ ! -e $CTRE_ESTI_ROOTED ]; then
    /opt/local/stow/Python3-3.8.1/bin/python3 $ROOTTREE \
        -i $CTRE_ESTI \
        -r "root" \
        -o  $CTRE_ESTI_ROOTED
    MYERRR=$(/opt/local/stow/Python3-3.8.1/bin/python3 $COMPARE -t1 $CTRE_ESTI -t2 $CTRE_ESTI_ROOTED)
    echo "$MYMODL,$MYMTHD,$MYERRR" > ${MYMTHD}x_debug_cell_lineage_tree_error.csv
fi

