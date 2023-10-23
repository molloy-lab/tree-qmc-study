#!/bin/bash

SCAL=$1
NGEN=$2
NBPS=$3
REPL=$4

MODL="$SCAL.$NGEN.$NBPS"
MYMODL="$SCAL,$NGEN,$NBPS,$REPL"

# Define directories
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJECTDIR="$GROUPDIR/ekmolloy/tree-qmc-study"
COMPARE="$PROJECTDIR/tools/compare_two_trees.py"

# Define input files
INDIR="$GROUPDIR/group/data/mahbub2021wqfm/37-taxon-mammalian-simulated/$MODL/$REPL"
STRE_TRUE="$INDIR/../../true_tree_trimmed"

# Define output files
OUTDIR="$PROJECTDIR/data/mahbub2021wqfm/37-taxon-mammalian-simulated/$MODL/$REPL"
GTRE_FILE="$OUTDIR/all_gt.tre"

if [ ! -e $INDIR/all_gt.tre ]; then
    echo "$INDIR/all_gt.tre does not exist!"
    exit
fi

if [ ! -e $GTRE_FILE ]; then
    cp $INDIR/all_gt.tre $GTRE_FILE
fi

MYMTHD="astral_3_v5.7.7"
MYSTRE="$OUTDIR/${MYMTHD}"
if [ ! -e $MYSTRE.tre ]; then
    ASTRAL3DIR="$GROUPDIR/group/software/ASTRAL_v5.7.7/Astral"
    ASTRAL3=astral.5.7.7.jar

    MYTIME="$(time (java -Xmx36G \
              -D"java.library.path=$ASTRAL3DIR/lib" \
              -jar $ASTRAL3DIR/$ASTRAL3 \
              -t0 \
              -i $GTRE_FILE \
              -o $MYSTRE.tre \
              &> $MYSTRE.log) 2>&1 1>/dev/null)"

    uname -a > ${MYSTRE}_node_info.csv
    MYNODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

    MYSECS="$(echo $MYTIME | awk '{print $2","$4","$6}')"
    echo "$MYMODL,$MYMTHD,$MYNODE,$MYSECS" > ${MYSTRE}_runtime.csv

    MYERRR=$(/opt/local/stow/Python3-3.8.1/bin/python3 $COMPARE -t1 $STRE_TRUE -t2 $MYSTRE.tre)
    echo "$MYMODL,$MYMTHD,$MYERRR" > ${MYSTRE}_species_tree_error.csv
fi

