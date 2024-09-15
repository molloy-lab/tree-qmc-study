#!/bin/bash

exit

# Define directories and files
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
DATADIR="$PROJDIR/data/morel2022asteroid-plants"
GTRE_FILE="all-relabeled-raxml-ng.LG+G.geneTree.newick"
WTREEQMC="$PROJDIR/software/TREE-QMC/tree-qmc"

cd $DATADIR

if [ ! -e $GTRE_FILE ]; then
    echo "ERROR: $GTRE_FILE does not exist!"
    exit 
fi

# Root species tree at Red algae for all methods
TREES=( "reference-speciesTree.newick" \
        "concatenation-single.LG+G.speciesTree.newick" \
        "treeqmc_wf_n2.tre" \
        "asteroid.bestTree.newick" \
	"wastrid_vanilla.tre" \
	"aster_v1.16.3.4.tre" )

END=${#TREES[@]}
END=$[END-1]

for INDEX in `seq 0 $END`; do
    MYSTRE=${TREES[$INDEX]}
    echo $MYSTRE
    echo

    if [ $MYSTRE == "wastrid_vanilla.tre" ]; then
        ROOT="ChondrusUUUCrispus"
    else
        ROOT="ChondrusUUUCrispus,CyanidioschyzonUUUMerolae,GaldieriaUUUSulphuraria"
    fi

    if [ -e $MYSTRE ] && [ ! -e rooted_${MYSTRE} ]; then
        $WTREEQMC \
            --root $ROOT \
            --rootonly $MYSTRE \
            -i $GTRE_FILE \
	    -o rooted_${MYSTRE}
    fi
done

