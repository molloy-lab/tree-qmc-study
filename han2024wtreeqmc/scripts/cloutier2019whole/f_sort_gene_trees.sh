#!/bin/bash

exit

TYPES=( "CNEE" "intron" "UCE" )

GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
DATADIR="$PROJDIR/data/cloutier2019whole"

for TYPE in ${TYPES[@]}; do
    cd $DATADIR/alignments/${TYPE}s
    echo "TYPE,GENE,NTAX,ALEN,NPIP,NTAX_nogallus,ALEN_nogallus,NPIP_nogallus" > tmp1.txt
    cat *csv > tmp2.txt
    cat tmp1.txt tmp2.txt ${TYPE}_alignment_info.csv
done


for TYPE in ${TYPES[@]}; do
    cd $DATADIR/RAxML_bestML_gene_trees/${TYPE}s
    pwd

    cat *.csv > tmp2.txt

    if [ ! -e ${TYPE}_branch_support.csv ]; then
        cat *.csv > ${TYPE}_branch_support.csv
    fi

    if [ ! -e ${TYPE}_branch_support_sorted.txt ]; then
        echo "TYPE,GENE,SAVG,SMED,SSTD,SMIN,SMAX" > tmp1.txt
        
        sed 's/,/ /g' tmp2.txt | \
            awk '{print $3,$2}' | \
            sort -r | \
            awk '{print NR-1,$2,$1}' \
            >> ${TYPE}_branch_support_sorted.txt

        cat tmp1.txt tmp2.txt > ${TYPE}_branch_support.csv
    fi

    rm tmp*.txt

    if [ ! -e ${TYPE}_abayes_gene_trees_sorted.tre ]; then
        GTRES=( $(awk '{print $2".tre.abayes"}' ${TYPE}_branch_support_sorted.txt) )
        cat ${GTRES[@]} > ${TYPE}_abayes_gene_trees_sorted.tre
    fi
done

