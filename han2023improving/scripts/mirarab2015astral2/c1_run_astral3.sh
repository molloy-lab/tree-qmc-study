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

MYMTHD="astral_3_v5.7.7"
MYSTRE="$OUTDIR/${MYMTHD}_${GTRE}_${NGEN}"
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

