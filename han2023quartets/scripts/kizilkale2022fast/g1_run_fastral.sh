#!/bin/bash

# Define software and utilities
GROUPDIR="/fs/cbcb-lab/ekmolloy"

PROJECTDIR="$GROUPDIR/ekmolloy/species-to-tumor-study"
COMPARE="$PROJECTDIR/tools/compare_two_trees.py"

MYMODL=$1
DATDIR=$2

#MYMTHD="fastral"
MYMTHD="fastral_wroot"

#INPUT="$DATDIR/noisy.nwk"
INPUT="$DATDIR/noisy_wroot.nwk"
CTRE_TRUE="$DATDIR/ground.cell_lineage_tree.nwk"
CTRE_BASE="$DATDIR/$MYMTHD"
CTRE_ESTI="$CTRE_BASE.cell_lineage_tree.nwk"
NTHREAD=1

if [ -e $CTRE_ESTI ]; then
    X=$(grep ";" $CTRE_ESTI | wc -l)
    if [ $X -eq 0 ]; then
        rm $CTRE_BASE*
    fi
fi

if [ ! -e $CTRE_ESTI ]; then
    ASTRAL3DIR="/cbcbhomes/ekmolloy/.local/lib/python3.8/site-packages/FASTRAL/ASTRAL-modified"
    ASTRAL3="astral.5.7.3.modified.jar"
    ASTRID2="/cbcbhomes/ekmolloy/.local/lib/python3.8/site-packages/FASTRAL/ASTRID/ASTRID-linux"
    FASTRAL="/cbcbhomes/ekmolloy/.local/bin/fastral"

    NGEN1=$(wc -l $INPUT | awk '{print $1}')
    NGEN2=$(echo "($NGEN1 / 2) + 1" | bc)
    NGEN4=$(echo "($NGEN1 / 4) + 1" | bc)
    NGEN8=$(echo "($NGEN1 / 8) + 1" | bc)

    MYTIME="$(time (python3 $FASTRAL \
              --ns 1,10,20,20 \
              --nt $NGEN1,$NGEN2,$NGEN4,$NGEN8 \
              --k $NGEN1 \
              --it $INPUT \
              --os "${CTRE_BASE}_output/samples" \
              --aggregate "${CTRE_BASE}_output/${MYMTHD}_constraints.tre" \
              --o "${CTRE_BASE}_output/${MYMTHD}.tre" \
              --time "${CTRE_BASE}_output/${MYMTHD}_time.log" \
              --path_ASTRAL $ASTRAL3DIR/$ASTRAL3 \
              --path_ASTRID $ASTRID2 \
              &> $CTRE_BASE.log) 2>&1 1>/dev/null)"

    cp ${CTRE_BASE}_output/${MYMTHD}.tre $CTRE_ESTI
    cp ${CTRE_BASE}_output/${MYMTHD}_constraints.tre ${CTRE_BASE}_constraints.tre
    cp ${CTRE_BASE}_output/${MYMTHD}_time.log ${CTRE_BASE}_time.log
    rm -rf ${CTRE_BASE}_output/

    uname -a > ${CTRE_BASE}_node_info.csv
    MYNODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

    MYSECS="$(echo $MYTIME | awk '{print $2","$4","$6}')"
    echo "$MYMODL,$MYMTHD,$MYNODE,$NTHREAD,$MYSECS" > ${CTRE_BASE}_runtime.csv

    MYERRR=$(python3 $COMPARE -t1 $CTRE_TRUE -t2 $CTRE_ESTI)
    echo "$MYMODL,$MYMTHD,$MYERRR" > ${CTRE_BASE}_cell_lineage_tree_error.csv
fi
