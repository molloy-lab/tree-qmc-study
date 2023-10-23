#!/bin/bash

# Define software and utilities
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJECTDIR="$GROUPDIR/ekmolloy/species-to-tumor-study"
TREEQMCBP="$PROJECTDIR/software/TREE-QMC-bip/TREE-QMC-bip"
ROOTTREE="$PROJECTDIR/tools/root_tree.py"
COMPARE="$PROJECTDIR/tools/compare_two_trees.py"

MYMODL=$1
DATDIR=$2

#INPUT="noisy.CFMatrix"       # Run 1
INPUT="noisy_wroot.CFMatrix"  # Run 2
CTRE_TRUE="ground.cell_lineage_tree.nwk"
NTHREAD=1

cd $DATDIR

#NS=( 0 1 2 ) # Run 1
NS=( 2 )      # Run 2

for N in ${NS}; do

MYMTHD="treeqmcbip_v1.0.0_n${N}_wroot"
CTRE_ESTI="$MYMTHD.cell_lineage_tree.nwk"

if [ -e $CTRE_ESTI ]; then
    X=$(grep ";" $CTRE_ESTI | wc -l)
    if [ $X -eq 0 ]; then
        rm ${MYMTHD}*
    fi
fi

if [ ! -e $CTRE_ESTI ]; then
    MYTIME="$(time ($TREEQMCBP -n $N \
                        -i $INPUT \
                        -o $CTRE_ESTI \
                        &> $MYMTHD.log) 2>&1 1>/dev/null)"

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

done

