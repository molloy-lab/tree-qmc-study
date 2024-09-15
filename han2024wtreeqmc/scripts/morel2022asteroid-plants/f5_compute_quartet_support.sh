#!/bin/bash

exit

# Define directories and files
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
ASTRAL3DIR="$GROUPDIR/group/software/ASTRAL_v5.7.7/Astral"
ASTRAL3=astral.5.7.7.jar
DATADIR="$PROJDIR/data/morel2022asteroid-plants"
GTRE_FILE="all-relabeled-raxml-ng.LG+G.geneTree.newick"
MYINFO="plants"

cd $DATADIR

if [ ! -e $GTRE_FILE ]; then
    echo "ERROR: $GTRE_FILE does not exist!"
    exit 
fi

# Score estimated species tree for all methods
MTHDS=( "reference" \
        "concatenation" \
        "treeqmc_wf_n2" \
        "asteroid" \
	"wastrid_vanilla" \
        "aster_v1.16.3.4" )

TREES=( "reference-speciesTree.newick" \
        "concatenation-single.LG+G.speciesTree.newick" \
        "treeqmc_wf_n2.tre" \
        "asteroid.bestTree.newick" \
	"wastrid_vanilla.tre" \
	"aster_v1.16.3.4.tre" )

END=${#TREES[@]}
END=$[END-1]

for INDEX in `seq 0 $END`; do
    MYMTHD=${MTHDS[$INDEX]}
    MYSTRE=${TREES[$INDEX]}
    echo $MYMTHD
    echo $MYSTRE
    echo

    if [ -e $MYSTRE ] && [ ! -e scored_${MYMTHD}.tre ]; then
        java -Xmx36G \
             -D"java.library.path=$ASTRAL3DIR/lib" \
             -jar $ASTRAL3DIR/$ASTRAL3 \
             -q $MYSTRE \
             -t2 \
             -i $GTRE_FILE \
             -o scored_${MYMTHD}.tre \
             &> scored_${MYMTHD}.log

        MYQSCR="$(grep "Final quartet score is:" scored_${MYMTHD}.log | awk '{print $5}')"
        MYNORM=$(grep "Final normalized quartet score" scored_${MYMTHD}.log | awk '{print $6}')
        echo "$MYINFO,$MYMTHD,$MYQSCR,$MYNORM" > ${MYMTHD}_quartet_score.csv
    fi
done


MTHDS=( "reference" )
TREES=( "reference-speciesTree.newick" )

END=${#TREES[@]}
END=$[END-1]

for INDEX in `seq 0 $END`; do
    MYMTHD=${MTHDS[$INDEX]}
    MYSTRE=${TREES[$INDEX]}
    echo $MYMTHD
    echo $MYSTRE
    echo

    if [ -e $MYSTRE ] && [ ! -e scored_t1_${MYMTHD}.tre ]; then
        java -Xmx36G \
             -D"java.library.path=$ASTRAL3DIR/lib" \
             -jar $ASTRAL3DIR/$ASTRAL3 \
             -q $MYSTRE \
             -t1 \
             -i $GTRE_FILE \
             -o scored_t1_${MYMTHD}.tre \
             &> scored_t1_${MYMTHD}.log

        MYQSCR="$(grep "Final quartet score is:" scored_t1_${MYMTHD}.log | awk '{print $5}')"
        MYNORM=$(grep "Final normalized quartet score" scored_t1_${MYMTHD}.log | awk '{print $6}')
        echo "$MYINFO,$MYMTHD,$MYQSCR,$MYNORM" > ${MYMTHD}_quartet_score.csv
    fi
done

