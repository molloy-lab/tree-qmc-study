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

MYMTHD="quartets"
MYQRTT="$OUTDIR/${MYMTHD}_${GTRE}_${NGEN}"
if [ ! -e $MYQRTT.qrt ]; then
    GET_QRTTS="$GROUPDIR/group/software/wQFM-2020/quartet-controller.sh"
    MYTIME="$(time (bash $GET_QRTTS \
        $GTRE_FILE \
        $MYQRTT.qrt \
        &> $MYQRTT.log) 2>&1 1>/dev/null)"

    uname -a > ${MYQRTT}_node_info.csv
    MYNODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

    MYSECS="$(echo $MYTIME | awk '{print $2","$4","$6}')"
    echo "$MYMODL,$MYMTHD,$MYNODE,$MYSECS" > ${MYQRTT}_runtime.csv
fi

QRTT_WQFM="$OUTDIR/wqfm_quartets_${GTRE}_${NGEN}"
if [ ! -e $QRTT_WQFM.qrt ]; then
    grep "),(" $MYQRTT.qrt >> $QRTT_WQFM.qrt
fi

QRTT_WQMC="$OUTDIR/wqmc_quartets_${GTRE}_${NGEN}"
if [ ! -e $QRTT_WQMC.qrt ]; then
    sed 's/),(/|/g' $QRTT_WQFM.qrt | \
        sed 's/(//g' | \
        sed 's/)//g' | \
        sed 's/; /:/g' \
        > $QRTT_WQMC.qrt
fi

