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
ASTRAL3DIR="$GROUPDIR/group/software/ASTRAL_v5.7.7/Astral"
ASTRAL3=astral.5.7.7.jar

# Define input files
INDIR="$GROUPDIR/group/data/mirarab2015astral2/$MODL/$REPL"
STRE_TRUE="$INDIR/s_tree.trees"

# Define output files
OUTDIR="$PROJECTDIR/data/$MODL/$REPL"
GTRE_FILE="$OUTDIR/${GTRE}_${NGEN}.txt"

if [ ! -e $GTRE_FILE ]; then
    head -n${NGEN} $INDIR/$GTRE > $GTRE_FILE
fi

echo "NTAX,PSIZE,SRATE,REPL,GTRE,NGEN,MTHD,QS,NQS" > $OUTDIR/header_quartet_score.csv

# Score the true species tree
MYSTRE="$OUTDIR/true_s_tree_${GTRE}_${NGEN}"
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
    MYSTRE="$OUTDIR/${MYMTHD}_${GTRE}_${NGEN}"
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

