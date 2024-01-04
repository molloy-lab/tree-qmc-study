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
ASTEROID="$GROUPDIR/group/software-compiled-on-EPYC-7313/Asteroid-nompi/Asteroid/build/bin/asteroid"

# Define input files
INDIR="$GROUPDIR/group/data/zhang2022weighting-dryad/S100/$REPL"
OUTDIR="$PROJDIR/data/zhang2022weighting-gteesim/S100/$REPL"
STRE_TRUE="$INDIR/s_tree.trees"
GTRE="$INDIR/bestMLestimatedgenetree/estimatedgenetre_${NBPS}.gtr.rerooted.final.contracted.non"
GTRE_FILE="estimatedgenetre.${NBPS}.${NGEN}"
if [ $SUPP == "abayes" ]; then
    GTRE="$GTRE.abayes"
    GTRE_FILE="estimatedgenetre.abayes.${NBPS}.${NGEN}"
    exit
fi

# Do work
cd $OUTDIR

if [ ! -e $GTRE_FILE ]; then
    head -n${NGEN} $GTRE > $GTRE_FILE
fi

MYMTHD="asteroid"
MYSTRE="${MYMTHD}_${SUPP}_${NBPS}bps_${NGEN}gen"
if [ ! -e $MYSTRE.bestTree.newick ]; then
    MYTIME="$(time ($ASTEROID \
                             -i $GTRE_FILE \
                             -p $MYSTRE \
                             &> $MYSTRE.log) 2>&1 1>/dev/null)"

    uname -a > ${MYSTRE}_node_info.csv
    MYNODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

    MYSECS="$(echo $MYTIME | awk '{print $2","$4","$6}')"
    echo "$MYMODL,$MYMTHD,$MYNODE,$MYSECS" > ${MYSTRE}_runtime.csv

    MYERRR=$(python3 $COMPARE -t1 $STRE_TRUE -t2 $MYSTRE.bestTree.newick)
    echo "$MYMODL,$MYMTHD,$MYERRR" > ${MYSTRE}_species_tree_error.csv
fi

