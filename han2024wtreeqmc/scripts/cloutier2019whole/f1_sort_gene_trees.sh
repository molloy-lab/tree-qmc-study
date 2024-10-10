#!/bin/bash

exit

TYPES=( "CNEEs" )
TYPES=( "introns" )
TYPES=( "UCEs_minus_bad105" "bad105" )

GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
DATADIR="$PROJDIR/data/cloutier2019whole/RAxML_bestML_gene_trees"

for TYPE in ${TYPES[@]}; do
    if [ $TYPE == "UCEs_minus_bad105" ]; then
        cd $DATADIR/UCEs
    elif [ $TYPE == "bad105" ]; then
        cd $DATADIR/UCEs/exclude
    else
        cd $DATADIR/${TYPE}
        pwd
    fi

    if [ ! -e ${TYPE}_branch_support_sorted.txt ]; then
        rm ${TYPE}_branch_support.csv
        cat *.csv > tmp2.txt

        echo "TYPE,GENE,SAVG,SMED,SSTD,SMIN,SMAX" > tmp1.txt
        
        sed 's/,/ /g' tmp2.txt | \
            awk '{print $3,$2}' | \
            sort -r | \
            awk '{print NR-1,$2,$1}' \
            >> ${TYPE}_branch_support_sorted.txt

        cat tmp1.txt tmp2.txt > ${TYPE}_branch_support.csv
        rm tmp*.txt
    fi

    if [ $TYPE == "bad105" ]; then
        if [ ! -e ${TYPE}_abayes_gene_trees_sorted.tre ]; then
            GTRES=( $(awk '{print "exclude."$2".tre.abayes"}' ${TYPE}_branch_support_sorted.txt) )
            cat ${GTRES[@]} > ${TYPE}_abayes_gene_trees_sorted.tre
        fi
    else
        if [ ! -e ${TYPE}_abayes_gene_trees_sorted.tre ]; then
            GTRES=( $(awk '{print $2".tre.abayes"}' ${TYPE}_branch_support_sorted.txt) )
            cat ${GTRES[@]} > ${TYPE}_abayes_gene_trees_sorted.tre
        fi
    fi
done

