#!/bin/bash

MODL=$1
REPL=$2


# Define directories
LABDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$LABDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
DATADIR="$PROJDIR/data/mirarab2015astral2-extsim/$MODL/$REPL"
GTRE="estimatedgenetre"
NGEN=1000

cd $DATADIR

if [ ! -e $GTRE.abayes ]; then
    GTRES=()

    for GTAG in `seq -f "%04g" 1 $NGEN`; do
        if [ ! -e $GTRE.abayes-$GTAG ]; then
            echo "$DATADIR/$GTRE.abayes-$GTAG does not exist!"
            exit 1
        fi

	TMP=$(grep ";" $GTRE.abayes-$GTAG) 
        if [ -z $TMP ]; then
            echo "$DATADIR/$GTRE.abayes-$GTAG does not contain a newick string!"
	    exit 1
        fi

	GTRES=( ${GTRES[@]} $GTRE.abayes-$GTAG )
    done

    cat ${GTRES[@]} > $GTRE.abayes
fi

if [ -e $GTRE.abayes ]; then
    TMP=$(grep ";" $GTRE.abayes | wc -l)
    if [ $TMP -ne $NGEN ]; then
        echo "$DATADIR/$GTRE.abayes contains wrong number of lines!"
	exit 1
    fi
fi

