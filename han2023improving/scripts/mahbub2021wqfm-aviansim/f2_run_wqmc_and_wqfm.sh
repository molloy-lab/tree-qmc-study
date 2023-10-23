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
COMPARE="$PROJECTDIR/tools/compare_two_trees.py"
RELABEL="$PROJECTDIR/tools/relabel_tree_by_list.py"

# Define input files
INDIR="$GROUPDIR/group/data/mahbub2021wqfm/48-taxon-avian-simulated/$MODL/$REPL"
STRE_TRUE="$INDIR/../../true_tree_trimmed"


# Define output files
OUTDIR="$PROJECTDIR/data/mahbub2021wqfm/48-taxon-avian-simulated/$MODL/$REPL"
QRTT_WQMC="$OUTDIR/wqmc_quartets"
QRTT_WQFM="$OUTDIR/wqfm_quartets"


MYMTHD="wqmc_v3.0"
MYSTRE="$OUTDIR/${MYMTHD}"
if [ ! -e $MYSTRE.tre ] && [ -e $QRTT_WQMC.qrt ] ; then
    WQMC="$GROUPDIR/group/software/wQMC_v3.0/max-cut-tree"

    MYTIME="$(time ($WQMC \
        qrtt=$QRTT_WQMC.qrt \
        weights=on \
        otre=${MYSTRE}_need_to_relabel.tre \
        &> $MYSTRE.log) 2>&1 1>/dev/null)"

    /opt/local/stow/Python3-3.8.1/bin/python3 $RELABEL \
        -t ${MYSTRE}_need_to_relabel.tre\
        -x ${QRTT_WQMC}_taxon_map.txt \
        -o $MYSTRE.tre

    uname -a > ${MYSTRE}_node_info.csv
    MYNODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

    MYSECS="$(echo $MYTIME | awk '{print $2","$4","$6}')"
    echo "$MYMODL,$MYMTHD,$MYNODE,$MYSECS" > ${MYSTRE}_runtime.csv

    MYERRR=$(/opt/local/stow/Python3-3.8.1/bin/python3 $COMPARE -t1 $STRE_TRUE -t2 $MYSTRE.tre)
    echo "$MYMODL,$MYMTHD,$MYERRR" > ${MYSTRE}_species_tree_error.csv
fi

MYMTHD="wqfm_v1.3"
MYSTRE="$OUTDIR/${MYMTHD}"
if [ ! -e $MYSTRE.tre ] && [ -e $QRTT_WQFM.qrt ]; then
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


