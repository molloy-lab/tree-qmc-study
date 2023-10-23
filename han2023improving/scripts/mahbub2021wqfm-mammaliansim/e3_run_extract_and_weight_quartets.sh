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
PREPARE="$PROJECTDIR/tools/prepare_quartets_for_wqmc.py"

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

MYMTHD="quartets"
MYQRTT="$OUTDIR/${MYMTHD}"
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

QRTT_WQFM="$OUTDIR/wqfm_quartets"
if [ ! -e $QRTT_WQFM.qrt ]; then
    grep "),(" $MYQRTT.qrt >> $QRTT_WQFM.qrt
fi

QRTT_WQMC="$OUTDIR/wqmc_quartets"
if [ ! -e $QRTT_WQMC.qrt ]; then
    /opt/local/stow/Python3-3.8.1/bin/python3 $PREPARE \
        -i $QRTT_WQFM.qrt \
        -o $QRTT_WQMC
fi

