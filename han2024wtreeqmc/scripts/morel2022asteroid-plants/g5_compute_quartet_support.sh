#!/bin/bash

#exit

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

TREES=( "rooted_reference-speciesTree.newick" \
        "rooted_concatenation-single.LG+G.speciesTree.newick" \
        "rooted_treeqmc_wf_n2.tre" \
        "rooted_asteroid.bestTree.newick" \
	"rooted_wastrid_vanilla.tre" \
	"rooted_aster_v1.16.3.4.tre" )

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
             -o scored_rooted_${MYMTHD}.tre \
             &> scored_rooted_${MYMTHD}.log

        MYQSCR="$(grep "Final quartet score is:" scored_rooted_${MYMTHD}.log | awk '{print $5}')"
        MYNORM=$(grep "Final normalized quartet score" scored_rooted_${MYMTHD}.log | awk '{print $6}')
        echo "$MYINFO,$MYMTHD,$MYQSCR,$MYNORM" > ${MYMTHD}_quartet_score.csv
    fi
done


MTHDS=( "reference" )
TREES=( "rooted_reference-speciesTree.newick" )

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
             -o scored_t1_rooted_${MYMTHD}.tre \
             &> scored_t1_rooted_${MYMTHD}.log

        MYQSCR="$(grep "Final quartet score is:" scored_t1_rooted_${MYMTHD}.log | awk '{print $5}')"
        MYNORM=$(grep "Final normalized quartet score" scored_t1_rooted_${MYMTHD}.log | awk '{print $6}')
        echo "$MYINFO,$MYMTHD,$MYQSCR,$MYNORM" > ${MYMTHD}_quartet_score.csv
    fi
done

#scored_rooted_asteroid.log:Final quartet score is: 1863661
#scored_rooted_aster_v1.16.3.4.log:Final quartet score is: 1868476
#scored_rooted_concatenation.log:Final quartet score is: 1861015
#scored_rooted_reference.log:Final quartet score is: 1717641
#scored_rooted_treeqmc_wf_n2.log:Final quartet score is: 1864686
#scored_rooted_wastrid_vanilla.log:Final quartet score is: 1842335
#scored_t1_rooted_reference.log:Final quartet score is: 1717641
