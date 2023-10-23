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

QRTT_WQMC="$OUTDIR/wqmc_quartets_${GTRE}_${NGEN}"
if [ ! -e $QRTT_WQMC ]; then
   echo "$QRTT_WQMC does not exist!"
fi

MYMTHD="wqmc_v3.0"
MYSTRE="$OUTDIR/${MYMTHD}_${GTRE}_${NGEN}"
if [ ! -e $MYSTRE.tre ]; then
    WQMC="$GROUPDIR/group/software/wQMC_v3.0/max-cut-tree"

    MYTIME="$(time ($WQMC \
        qrtt=$QRTT_WQMC.qrt \
        weights=on \
        otre=$MYSTRE.tre \
        &> $MYSTRE.log) 2>&1 1>/dev/null)"

    uname -a > ${MYSTRE}_node_info.csv
    MYNODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

    MYSECS="$(echo $MYTIME | awk '{print $2","$4","$6}')"
    echo "$MYMODL,$MYMTHD,$MYNODE,$MYSECS" > ${MYSTRE}_runtime.csv

    MYERRR=$(/opt/local/stow/Python3-3.8.1/bin/python3 $COMPARE -t1 $STRE_TRUE -t2 $MYSTRE.tre)
    echo "$MYMODL,$MYMTHD,$MYERRR" > ${MYSTRE}_species_tree_error.csv
fi

QRTT_WQFM="$OUTDIR/wqfm_quartets_${GTRE}_${NGEN}"
if [ ! -e $QRTT_WQFM ]; then
   echo "$QRTT_WQFM does not exist!"
fi

MYMTHD="wqfm_v1.3"
MYSTRE="$OUTDIR/${MYMTHD}_${GTRE}_${NGEN}"
if [ ! -e $MYSTRE.tre ]; then
    WQFM="$GROUPDIR/group/software/wQFM-2020/wQFM-v1.3.jar"

    MYTIME="$(time (java -Xmx36G \
        -jar $WQFM \
        -i $QRTT_WQFM.qrt \
        -o $MYSTRE.tre \
        &> $MYSTRE.log) 2>&1 1>/dev/null)"

    uname -a > ${MYSTRE}_node_info.csv
    MYNODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

    MYSECS="$(echo $MYTIME | awk '{print $2","$4","$6}')"
    echo "$MYMODL,$MYMTHD,$MYNODE,$MYSECS" > ${MYSTRE}_runtime.csv

    MYERRR=$(/opt/local/stow/Python3-3.8.1/bin/python3 $COMPARE -t1 $STRE_TRUE -t2 $MYSTRE.tre)
    echo "$MYMODL,$MYMTHD,$MYERRR" > ${MYSTRE}_species_tree_error.csv
fi


