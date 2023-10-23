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
TREEQMC="$PROJECTDIR/software/TREE-QMC/TREE-QMC"

# Define input files
INDIR="$GROUPDIR/group/data/mahbub2021wqfm/48-taxon-avian-simulated/$MODL/$REPL"
STRE_TRUE="$INDIR/../../true_tree_trimmed"

# Define output files
OUTDIR="$PROJECTDIR/data/mahbub2021wqfm/48-taxon-avian-simulated/$MODL/$REPL"
GTRE_FILE="$OUTDIR/all_gt.tre"

if [ ! -e $INDIR/all_gt.tre ]; then
    echo "$INDIR/all_gt.tre does not exist!"
    exit
fi

if [ ! -e $GTRE_FILE ]; then
    cp $INDIR/all_gt.tre $GTRE_FILE
fi

MYTIME="$(time (ls) 2>&1 1>/dev/null)"
MYDATA="$(echo $MYTIME | awk '{print $1","$3","$5}')"
echo "SCAL,NGEN,NBPS,REPL,MTHD,NODE,$MYDATA" > $OUTDIR/header_runtime.csv
echo "SCAL,NGEN,NBPS,REPL,MTHD,NL,I1,I2,FN,FP,RF" > $OUTDIR/header_species_tree_error.csv

# Run all methods

MYMTHD="treeqmc_n0_v1.0.0"
MYSTRE="$OUTDIR/${MYMTHD}"
if [ ! -e $MYSTRE.tre ]; then
    MYTIME="$(time ($TREEQMC -n 0 \
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

MYMTHD="treeqmc_n1_v1.0.0"
MYSTRE="$OUTDIR/${MYMTHD}"
if [ ! -e $MYSTRE.tre ]; then
    if [ -e $GTRE_FILE.refined ]; then
    MYTIME="$(time ($TREEQMC -n 1 \
                             -i $GTRE_FILE.refined \
                             -o $MYSTRE.tre \
                             &> $MYSTRE.log) 2>&1 1>/dev/null)"
    else
    MYTIME="$(time ($TREEQMC -n 1 \
                             -i $GTRE_FILE \
                             -o $MYSTRE.tre \
                             &> $MYSTRE.log) 2>&1 1>/dev/null)"
    fi

    uname -a > ${MYSTRE}_node_info.csv
    MYNODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

    MYSECS="$(echo $MYTIME | awk '{print $2","$4","$6}')"
    echo "$MYMODL,$MYMTHD,$MYNODE,$MYSECS" > ${MYSTRE}_runtime.csv

    MYERRR=$(/opt/local/stow/Python3-3.8.1/bin/python3 $COMPARE -t1 $STRE_TRUE -t2 $MYSTRE.tre)
    echo "$MYMODL,$MYMTHD,$MYERRR" > ${MYSTRE}_species_tree_error.csv
fi

MYMTHD="treeqmc_n2_v1.0.0"
MYSTRE="$OUTDIR/${MYMTHD}"
if [ ! -e $MYSTRE.tre ]; then
    if [ -e $GTRE_FILE.refined ]; then
    MYTIME="$(time ($TREEQMC -n 2 \
                             -i $GTRE_FILE.refined \
                             -o $MYSTRE.tre \
                             &> $MYSTRE.log) 2>&1 1>/dev/null)"
    else
    MYTIME="$(time ($TREEQMC -n 2 \
                             -i $GTRE_FILE \
                             -o $MYSTRE.tre \
                             &> $MYSTRE.log) 2>&1 1>/dev/null)"
    fi

    uname -a > ${MYSTRE}_node_info.csv
    MYNODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

    MYSECS="$(echo $MYTIME | awk '{print $2","$4","$6}')"
    echo "$MYMODL,$MYMTHD,$MYNODE,$MYSECS" > ${MYSTRE}_runtime.csv

    MYERRR=$(/opt/local/stow/Python3-3.8.1/bin/python3 $COMPARE -t1 $STRE_TRUE -t2 $MYSTRE.tre)
    echo "$MYMODL,$MYMTHD,$MYERRR" > ${MYSTRE}_species_tree_error.csv
fi

