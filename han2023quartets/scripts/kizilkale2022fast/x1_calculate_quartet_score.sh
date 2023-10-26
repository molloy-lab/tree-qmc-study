#!/bin/bash

GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJECTDIR="$GROUPDIR/ekmolloy/species-to-tumor-study"

ASTRAL3DIR="$GROUPDIR/group/software/ASTRAL_v5.7.7/Astral"
ASTRAL3=astral.5.7.7.jar

MYMODL=$1
DATDIR=$2

cd $DATDIR

pwd

MTHDS=( "huntress_v0.1.2.0_default" \
        "treeqmcbip_v1.0.0_n2_wrootx" \
        "scistree_v1.2.0.6" \
        "fastme_v2.1.5_wrootx" \
        "fastral_wrootx" )



for MYMTHD in ${MTHDS[@]}; do
    CTRE_ESTI="${MYMTHD}.cell_lineage_tree"
    DATA=""
    for INPUT in "noisy" "ground"; do
        if [ -e $CTRE_ESTI.nwk ] && [ ! -e ${CTRE_ESTI}_scored_${INPUT}.nwk ]; then
            java -Xmx36G \
                -D"java.library.path=$ASTRAL3DIR/lib" \
                -jar $ASTRAL3DIR/$ASTRAL3 \
                -q $CTRE_ESTI.nwk \
                -t1 \
                -i $INPUT.nwk \
                -o ${CTRE_ESTI}_scored_${INPUT}.nwk \
                &> ${CTRE_ESTI}_scored_${INPUT}.log
        fi

        NOISY_MYQSCR="$(grep "Final quartet score is:" ${CTRE_ESTI}_scored_${INPUT}.log | awk '{print $5}')"
        NOISY_MYNORM=$(grep "Final normalized quartet score" ${CTRE_ESTI}_scored_${INPUT}.log | awk '{print $6}')

        DATA="$DATA,$NOISY_MYQSCR,$NOISY_MYNORM"
    done
    echo "$MYMODL,${MYMTHD}${DATA}" > ${CTRE_ESTI}_quartet_score.csv
done

