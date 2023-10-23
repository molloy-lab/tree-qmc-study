#!/bin/bash

NTAX=$1
PSIZE=$2
SRATE=$3
REPL=$4
GTRE=$5
NGEN=$6

MODL="model.$NTAX.$PSIZE.$SRATE"
MYMODL="$NTAX,$PSIZE,$SRATE,$REPL,$GTRE,$NGEN"

# Define directories
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJECTDIR="$GROUPDIR/ekmolloy/tree-qmc-study"
COMPARE="$PROJECTDIR/tools/compare_two_trees.py"

# Define input files
INDIR="$GROUPDIR/group/data/mirarab2015astral2/$MODL/$REPL"
STRE_TRUE="$INDIR/s_tree.trees"

# Define output files
OUTDIR="$PROJECTDIR/data/mirarab2015astral2/$MODL/$REPL"
GTRE_FILE="$OUTDIR/${GTRE}_${NGEN}.txt"

if [ ! -e $GTRE_FILE ]; then
    head -n${NGEN} $INDIR/$GTRE > $GTRE_FILE
fi

MYMTHD="fastral"
MYSTRE="$OUTDIR/${MYMTHD}_${GTRE}_${NGEN}"
if [ ! -e $MYSTRE.tre ]; then
    ASTRAL3DIR="/cbcbhomes/ekmolloy/.local/lib/python3.8/site-packages/FASTRAL/ASTRAL-modified"
    ASTRAL3="astral.5.7.3.modified.jar"
    ASTRID2="/cbcbhomes/ekmolloy/.local/lib/python3.8/site-packages/FASTRAL/ASTRID/ASTRID-linux"
    FASTRAL="/cbcbhomes/ekmolloy/.local/bin/fastral"

    if [ $NGEN -eq 1000 ]; then
        FANT="1000,500,250,100"
    elif [ $NGEN -eq 250 ]; then
        FANT="250,125,63,25"
    else
        echo "Unable to set up FASTRAL!"
        exit
    fi
    echo "Running FASTRAL with $NGEN and $FANT"

    MYTIME="$(time ($FASTRAL \
              --ns 1,10,20,20 \
              --nt $FANT \
              --k $NGEN \
              --it $GTRE_FILE \
              --os "${MYSTRE}_output/samples" \
              --aggregate "${MYSTRE}_output/${MYMTHD}_${GTRE}_${NGEN}_constraints.tre" \
              --o "${MYSTRE}_output/${MYMTHD}_${GTRE}_${NGEN}.tre" \
              --time "${MYSTRE}_output/${MYMTHD}_${GTRE}_${NGEN}_time.log" \
              --path_ASTRAL $ASTRAL3DIR/$ASTRAL3 \
              --path_ASTRID $ASTRID2 \
              &> $MYSTRE.log) 2>&1 1>/dev/null)"

    cp ${MYSTRE}_output/${MYMTHD}_${GTRE}_${NGEN}.tre $MYSTRE.tre
    cp ${MYSTRE}_output/${MYMTHD}_${GTRE}_${NGEN}_constraints.tre ${MYSTRE}_constraints.tre
    cp ${MYSTRE}_output/${MYMTHD}_${GTRE}_${NGEN}_time.log ${MYSTRE}_time.log
    rm -rf ${MYSTRE}_output/

    uname -a > ${MYSTRE}_node_info.csv
    MYNODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

    MYTIME="$(echo $MYTIME | awk '{print $2","$4","$6}')"
    echo "$MYMODL,$MYMTHD,$MYNODE,$MYTIME" > ${MYSTRE}_runtime.csv

    MYERRR=$(/opt/local/stow/Python3-3.8.1/bin/python3 $COMPARE -t1 $STRE_TRUE -t2 $MYSTRE.tre)
    echo "$MYMODL,$MYMTHD,$MYERRR" > ${MYSTRE}_species_tree_error.csv
fi

