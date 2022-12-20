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
ASTRAL3DIR="$GROUPDIR/group/software/ASTRAL_v5.7.7/Astral"
ASTRAL3=astral.5.7.7.jar

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

echo "SCAL,NGEN,NBPS,REPL,MTHD,QS,NQS" > $OUTDIR/header_quartet_score.csv

# Score the true species tree
MYSTRE="$OUTDIR/true_s_tree"
if [ ! -e ${MYSTRE}_scored.tre ]; then
    java -Xmx36G \
         -D"java.library.path=$ASTRAL3DIR/lib" \
         -jar $ASTRAL3DIR/$ASTRAL3 \
         -q $STRE_TRUE \
         -t2 \
         -i $GTRE_FILE \
         -o ${MYSTRE}_scored.tre \
         &> ${MYSTRE}_scored.log

    MYQSCR="$(grep "Final quartet score is:" ${MYSTRE}_scored.log | awk '{print $5}')"
    MYNORM=$(grep "Final normalized quartet score" ${MYSTRE}_scored.log | awk '{print $6}')

    echo "$MYMODL,TRUE,$MYQSCR,$MYNORM" > ${MYSTRE}_quartet_score.csv
fi

# Score estimated species tree for all methods
MTHDS=( "treeqmc_n0_v1.0.0" \
        "treeqmc_n1_v1.0.0" \
        "treeqmc_n2_v1.0.0" \
        "astral_3_v5.7.7" \
        "fastral" \
        "wqfm_v1.3" \
        "wqmc_v3.0" )

for MYMTHD in ${MTHDS[@]}; do
    MYSTRE="$OUTDIR/${MYMTHD}"
    if [ -e $MYSTRE.tre ] && [ ! -e ${MYSTRE}_scored.tre ]; then
        java -Xmx36G \
            -D"java.library.path=$ASTRAL3DIR/lib" \
            -jar $ASTRAL3DIR/$ASTRAL3 \
            -q $MYSTRE.tre \
            -t2 \
            -i $GTRE_FILE \
            -o ${MYSTRE}_scored.tre \
            &> ${MYSTRE}_scored.log

        MYQSCR="$(grep "Final quartet score is:" ${MYSTRE}_scored.log | awk '{print $5}')"
        MYNORM=$(grep "Final normalized quartet score" ${MYSTRE}_scored.log | awk '{print $6}')

        echo "$MYMODL,$MYMTHD,$MYQSCR,$MYNORM" > ${MYSTRE}_quartet_score.csv
    fi
done

