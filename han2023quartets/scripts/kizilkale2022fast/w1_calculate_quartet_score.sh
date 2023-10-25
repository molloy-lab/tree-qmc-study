#!/bin/bash

GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJECTDIR="$GROUPDIR/ekmolloy/species-to-tumor-study"

ASTRAL3DIR="$GROUPDIR/group/software/ASTRAL_v5.7.7/Astral"
ASTRAL3=astral.5.7.7.jar

MYMODL=$1
DATDIR=$2

cd $DATDIR

pwd

INPUT="$DATDIR/noisy.nwk"

MTHDS=( "huntress_v0.1.2.0_default" \
        "treeqmcbip_v1.0.0_n2_wrootx" \
        "scistree_v1.2.0.6" \
        "fastme_v2.1.5_wrootx" \
        "fastral_wrootx" )

for MYMTHD in ${MTHDS[@]}; do
    CTRE_ESTI="${MYMTHD}.cell_lineage_tree"

    MYSTRE="$OUTDIR/${MYMTHD}"
    if [ -e $CTRE_ESTI.nwk ] && [ ! -e ${CTRE_ESTI}_scored.nwk ]; then
        java -Xmx36G \
            -D"java.library.path=$ASTRAL3DIR/lib" \
            -jar $ASTRAL3DIR/$ASTRAL3 \
            -q $CTRE_ESTI.nwk \
            -t1 \
            -i $INPUT \
            -o ${CTRE_ESTI}_scored.nwk \
            &> ${CTRE_ESTI}_scored.log

        MYQSCR="$(grep "Final quartet score is:" ${CTRE_ESTI}_scored.log | awk '{print $5}')"
        MYNORM=$(grep "Final normalized quartet score" ${CTRE_ESTI}_scored.log | awk '{print $6}')

        echo "$MYMODL,$MYMTHD,$MYQSCR,$MYNORM" > ${CTRE_ESTI}_quartet_score.csv
    fi
done
