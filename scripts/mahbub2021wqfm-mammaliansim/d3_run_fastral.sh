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

MYMTHD="fastral"
MYSTRE="$OUTDIR/${MYMTHD}"
if [ ! -e $MYSTRE.tre ]; then
    ASTRAL3DIR="/cbcbhomes/ekmolloy/.local/lib/python3.8/site-packages/FASTRAL/ASTRAL-modified"
    ASTRAL3="astral.5.7.3.modified.jar"
    ASTRID2="/cbcbhomes/ekmolloy/.local/lib/python3.8/site-packages/FASTRAL/ASTRID/ASTRID-linux"
    FASTRAL="/cbcbhomes/ekmolloy/.local/bin/fastral"

    XNGEN=$(echo $NGEN | sed 's/g//g')
    if [ $XNGEN -eq 800 ]; then
        FANT="800,400,200,100"
    elif [ $XNGEN -eq 400 ]; then
        FANT="400,200,100,50"
    elif [ $XNGEN -eq 200 ]; then
        FANT="200,100,50,25"
    elif [ $XNGEN -eq 100 ]; then
        FANT="100,50,25,13"
    elif [ $XNGEN -eq 50 ]; then
        FANT="50,25,13,7"
    elif [ $XNGEN -eq 25 ]; then
        FANT="25,13,7,4"
    else
        echo "Unable to set up FASTRAL!"
        exit
    fi
    echo "Running FASTRAL with $XNGEN and $FANT"

    MYTIME="$(time ($FASTRAL \
              --ns 1,10,20,20 \
              --nt $FANT \
              --k $XNGEN \
              --it $GTRE_FILE \
              --os "${MYSTRE}_output/samples" \
              --aggregate "${MYSTRE}_output/${MYMTHD}_constraints.tre" \
              --o "${MYSTRE}_output/${MYMTHD}.tre" \
              --time "${MYSTRE}_output/${MYMTHD}_time.log" \
              --path_ASTRAL $ASTRAL3DIR/$ASTRAL3 \
              --path_ASTRID $ASTRID2 \
              &> $MYSTRE.log) 2>&1 1>/dev/null)"

    cp ${MYSTRE}_output/${MYMTHD}.tre $MYSTRE.tre
    cp ${MYSTRE}_output/${MYMTHD}_constraints.tre ${MYSTRE}_constraints.tre
    cp ${MYSTRE}_output/${MYMTHD}_time.log ${MYSTRE}_time.log
    rm -rf ${MYSTRE}_output/

    uname -a > ${MYSTRE}_node_info.csv
    MYNODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

    MYTIME="$(echo $MYTIME | awk '{print $2","$4","$6}')"
    echo "$MYMODL,$MYMTHD,$MYNODE,$MYTIME" > ${MYSTRE}_runtime.csv

    MYERRR=$(/opt/local/stow/Python3-3.8.1/bin/python3 $COMPARE -t1 $STRE_TRUE -t2 $MYSTRE.tre)
    echo "$MYMODL,$MYMTHD,$MYERRR" > ${MYSTRE}_species_tree_error.csv
fi

