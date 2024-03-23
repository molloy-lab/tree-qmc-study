#!/bin/bash

REPL=$1
NBPS=$2
SUPP=$3
NGEN=$4

MYMODL="$REPL,$NBPS,$SUPP,$NGEN"

# Define directories
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"

# Define software
COMPARE="$PROJDIR/tools/compare_two_trees.py"
WTREEQMC="$PROJDIR/software/TREE-QMC/tree-qmc"

# Define input files
INDIR="$GROUPDIR/group/data/zhang2022weighting-dryad/S100/$REPL"
OUTDIR="$PROJDIR/data/zhang2022weighting-gteesim/S100/$REPL"
STRE_TRUE="$INDIR/s_tree.trees"
GTRE="$INDIR/bestMLestimatedgenetree/estimatedgenetre_${NBPS}.gtr.rerooted.final.contracted.non"
GTRE_FILE="estimatedgenetre.${NBPS}.${NGEN}"
ROPTS="--bootstrap"
if [ $SUPP == "abayes" ]; then
    ROPTS="--bayes"
    GTRE="$GTRE.abayes"
    GTRE_FILE="estimatedgenetre.abayes.${NBPS}.${NGEN}"
fi

# Do work
cd $OUTDIR

if [ ! -e $GTRE_FILE ]; then
    head -n${NGEN} $GTRE > $GTRE_FILE
fi

MYMTHDS=( "wtreeqmc_ws_n2" \
	  "wtreeqmc_wh_n1" \
	  "wtreeqmc_wh_n0" )

for MYMTHD in ${MYMTHDS[@]}; do
    if [ $MYMTHD == "wtreeqmc_ws_n2" ]; then
        OPTS="-w s $ROPTS"
    elif [ $MYMTHD == "wtreeqmc_wh_n1" ]; then
        OPTS="-w h --norm_atax 1 $ROPTS"
    elif [ $MYMTHD == "wtreeqmc_wh_n0" ]; then
        OPTS="-w h --norm_atax 0 $ROPTS"
    else
        echo "Do not recognize $MYMTHD"
	exit
    fi

    MYSTRE="${MYMTHD}_${SUPP}_${NBPS}bps_${NGEN}gen"
    if [ ! -e $MYSTRE.tre ]; then
        MYTIME="$(time ($WTREEQMC $OPTS \
                             -i $GTRE_FILE \
                             -o $MYSTRE.tre \
                             &> $MYSTRE.log) 2>&1 1>/dev/null)"

        uname -a > ${MYSTRE}_node_info.csv
        MYNODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

        MYSECS="$(echo $MYTIME | awk '{print $2","$4","$6}')"
        echo "$MYMODL,$MYMTHD,$MYNODE,$MYSECS" > ${MYSTRE}_runtime.csv

        MYERRR=$(python3 $COMPARE -t1 $STRE_TRUE -t2 $MYSTRE.tre)
        echo "$MYMODL,$MYMTHD,$MYERRR" > ${MYSTRE}_species_tree_error.csv
    fi
done

