#!/bin/bash

#exit

# ILS STUDY PARAMETERS
#NTAXS=( 200 )                            # Number of taxa
#HGHTS=( "10000000" "2000000" "500000" )  # Species tree height (number of generations)
#RATES=( "0.0000001" "0.000001" )         # Speciation rate
#REPLS=( $(seq  -f "%02g" 1 50) )         # Replicates

# SCALABILITY STUDY PARAMETERS
NTAXS=( 10 50 100 500 1000 )      # Number of taxa
NTAXS=( 500 1000 )
HGHTS=( "2000000" )               # Species tree height (number of generations)
RATES=( "0.000001" )              # Speciation rate
REPLS=( $(seq  -f "%02g" 1 50) )  # Replicates

for NTAX in ${NTAXS[@]}; do
    for HGHT in ${HGHTS[@]}; do
        for RATE in ${RATES[@]}; do
            MODL="model.$NTAX.$HGHT.$RATE"
            for REPL in ${REPLS[@]}; do
                if [ $NTAX -eq 1000 ]; then
                    GSS=( 1 101 201 301 401 501 601 701 801 901 )
	            GES=( 100 200 300 400 500 600 700 800 900 1000 )
                elif [ $NTAX -eq 500 ]; then
                    GSS=( 1 201 401 601 801 )
	            GES=( 200 400 600 800 1000 )
                elif [ $NTAX -eq 200 ]; then
                    GSS=( 1 251 501 751 )
	            GES=( 250 500 750 100 )
                else
                    GSS=( 1 )
                    GES=( 1000 )
                fi
                NBATCH=${#GSS[@]}
                for BIND in `seq 0 $[NBATCH-1]`; do
                    GS=${GSS[$BIND]}
                    GE=${GES[$BIND]}
     
echo "Submitting $MODL $REPL $GS $GE..."
sbatch \
    --job-name="b1.$MODL.$REPL.$GS.$GE" \
    --output="b1.$MODL.$REPL.$GS.$GE.%j.out" \
    --error="b1.$MODL.$REPL.$GS.$GE.%j.err" \
    --export=NTAX="$NTAX",HGHT="$HGHT",RATE="$RATE",REPL="$REPL",GS="$GS",GE="$GE" \
b1_drive.sbatch

                done
            done
        done
    done
done

