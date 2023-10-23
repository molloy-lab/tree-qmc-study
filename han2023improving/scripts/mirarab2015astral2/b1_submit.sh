#!/bin/bash

exit

# ILS
#NTAXS=( 200 )                             # Number of taxa
#PSIZES=( "10000000" "2000000" "500000" )  # Population size
#SRATES=( "0.0000001" "0.000001" )         # Speciation rate
#REPLS=( $(seq  -f "%02g" 1 50) )          # Replicates
#GTRE="estimatedgenetre"                   # Use estimated gene trees
##GTRE="truegenetrees"          # Use true gene trees
#NGENS=( 250 1000 )                        # Number of gene trees

# SCALABILITY
#NTAXS=( 10 50 100 500 1000 )      # Number of taxa
#PSIZES=( "2000000" )              # Population size
#SRATES=( "0.000001" )             # Speciation rate
#REPLS=( $(seq  -f "%02g" 1 50) )  # Replicates
#GTRE="estimatedgenetre"           # Use estimated gene trees
##GTRE="truegenetrees"          # Use true gene trees
#NGENS=( 250 1000 )                # Number of gene trees

for NTAX in ${NTAXS[@]}; do
    for PSIZE in ${PSIZES[@]}; do
        for SRATE in ${SRATES[@]}; do
            MODL="model.$NTAX.$PSIZE.$SRATE"
            for NGEN in ${NGENS[@]}; do
                for REPL in ${REPLS[@]}; do
echo "Submitting $MODL $REPL..."
sbatch \
    --job-name="b1.$MODL.$REPL.$GTRE.$NGEN" \
    --output="b1.$MODL.$REPL.$GTRE.$NGEN.%j.out" \
    --error="b1.$MODL.$REPL.$GTRE.$NGEN.%j.err" \
    --export=NTAX="$NTAX",PSIZE="$PSIZE",SRATE="$SRATE",REPL="$REPL",GTRE="$GTRE",NGEN="$NGEN" \
b1_drive.sbatch
                done
            done
        done
    done
done

