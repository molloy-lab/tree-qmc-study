#!/bin/bash

SCAL=$1
NGEN=$2
NBPS=$3
REPL=$4

MODL="$SCAL-$NGEN-$NBPS"
MYMODL="$SCAL,$NGEN,$NBPS,$REPL"

# Define directories
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJECTDIR="$GROUPDIR/ekmolloy/tree-qmc-study"
WQFIT="$PROJECTDIR/tools/wQfit.py"

# Define input files
INDIR="$GROUPDIR/group/data/mahbub2021wqfm/48-taxon-avian-simulated/$MODL/$REPL"
STRE_TRUE="$INDIR/../../true_tree_trimmed"

# Define output files
OUTDIR="$PROJECTDIR/data/mahbub2021wqfm/48-taxon-avian-simulated/$MODL/$REPL"
QRTT_WQFM="$OUTDIR/wqfm_quartets.qrt"

if [ ! -e $QRTT_WQFM ]; then
    echo "$QRTT_WQFM does not exist!"
    exit
fi

cd $OUTDIR

# Score estimated species tree for all methods
MTHDS=( "treeqmc_n0_v1.0.0" \
        "treeqmc_n1_v1.0.0" \
        "treeqmc_n2_v1.0.0" \
        "astral_3_v5.7.7" \
        "fastral" \
        "wqfm_v1.3" \
        "wqmc_v3.0" )

for MYMTHD in ${MTHDS[@]}; do
    STRE_ESTI="$MYMTHD.tre"
    OUTF="${MYMTHD}_wqfit.csv"

    if [ -e $STRE_ESTI ] && [ ! -e $OUTF ]; then
        MYQFIT=$(/opt/local/stow/Python3-3.8.1/bin/python3 $WQFIT -t1 $STRE_TRUE -t2 $STRE_ESTI -qw $QRTT_WQFM)
        echo "$MYMODL,$MYMTHD,$MYQFIT" > $OUTF
    fi
done

