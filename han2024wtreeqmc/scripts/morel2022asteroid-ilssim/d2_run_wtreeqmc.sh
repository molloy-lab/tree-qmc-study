#!/bin/bash

# Define software and utilities
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"

REFINE="$PROJDIR/tools/randomly_refine_trees.py"
COMPARE="$PROJDIR/tools/compare_two_trees.py"
WTREEQMC="$PROJDIR/software/TREE-QMC/tree-qmc"

MYMODL=$1
OUTDIR=$2
GTRE_FILE=$3
STRE_TRUE=$4

cd $OUTDIR

if [ ! -e $GTRE_FILE ]; then
    echo "ERROR: $GTRE_FILE does not exist!"
    exit
fi

MYMTHDS=( "wtreeqmc_wf_n0" \
	  "wtreeqmc_wf_n1" \
	  "wtreeqmc_wf_n2" \
	  "wtreeqmc_wf_n1_shared" \
	  "wtreeqmc_wf_n2_shared" \
	  "wtreeqmc_wn_n2" )

for MYMTHD in ${MYMTHDS[@]}; do
    if [ $MYMTHD == "wtreeqmc_wf_n0" ]; then
        OPTS="-w f --norm_atax 0"
    elif [ $MYMTHD == "wtreeqmc_wf_n1" ]; then
        OPTS="-w f --norm_atax 1"
    elif [ $MYMTHD == "wtreeqmc_wf_n2" ]; then
        OPTS="-w f --norm_atax 2"
    elif [ $MYMTHD == "wtreeqmc_wf_n1_shared" ]; then
        OPTS="-w f --norm_atax 1 --shared"
    elif [ $MYMTHD == "wtreeqmc_wf_n2_shared" ]; then
        OPTS="-w f --norm_atax 2 --shared"
    elif [ $MYMTHD == "wtreeqmc_wn_n2" ]; then
        OPTS="-w n --norm_atax 2"
    else
        echo "Do not recognize $MYMTHD"
	exit
    fi

    if [ ! -e $MYMTHD.tre ]; then
        MYTIME="$(time ($WTREEQMC $OPTS -i $GTRE_FILE -o $MYMTHD.tre &> $MYMTHD.log) 2>&1 1>/dev/null)"

        uname -a > ${MYMTHD}_node_info.csv
        MYNODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

        MYSECS="$(echo $MYTIME | awk '{print $2","$4","$6}')"
        echo "$MYMODL,$MYMTHD,$MYNODE,$MYSECS" > ${MYMTHD}_runtime.csv

        MYERRR=$(python3 $COMPARE -t1 $STRE_TRUE -t2 $MYMTHD.tre)
        echo "$MYMODL,$MYMTHD,$MYERRR" > ${MYMTHD}_species_tree_error.csv
    fi

    if [ ! -e ${MYMTHD}_refined.tre ]; then
        python3 $REFINE -i $MYMTHD.tre &> ${MYMTHD}_refined.tre
    
        MYERRR=$(python3 $COMPARE -t1 $STRE_TRUE -t2 ${MYMTHD}_refined.tre)
        echo "$MYMODL,${MYMTHD}_refined,$MYERRR" > ${MYMTHD}_refined_species_tree_error.csv
    fi
done

