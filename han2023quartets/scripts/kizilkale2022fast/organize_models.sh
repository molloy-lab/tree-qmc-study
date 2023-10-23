#!/bin/bash

exit

NCELL=( "n_100" "n_200" "n_300" "n_1000" "n_5000" )
NMUTS=( "m_50" "m_100" "m_200" "m_300" "m_500" "m_1000" )

FILES=()
for NC in ${NCELL[@]}; do
    for NM in ${NMUTS[@]}; do
        grep "$NC-" model_list.txt | grep "$NM-" > model-list-${NC}-${NM}.txt
        FILES=( ${FILES[@]} "model-list-${NC}-${NM}.txt" )
    done
done

echo -n "" > model_list_sorted.txt
for X in ${FILES[@]}; do
    #echo $X
    cat $X >> model_list_sorted.txt
done

rm model-list-n_*

