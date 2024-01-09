#!/bin/bash

NTAX=$1
HGHT=$2
RATE=$3
REPL=$4
SUPP=$5
NGEN=$6

MODL="model.$NTAX.$HGHT.$RATE"
MYMODL="$NTAX,$HGHT,$RATE,$REPL,$SUPP,$NGEN"

# Define directories
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"

# Define software
COMPARE="$PROJDIR/tools/compare_two_trees.py"
ASTER="$GROUPDIR/group/software-compiled-on-EPYC-7313/ASTER/bin/astral"
ASTERH="$GROUPDIR/group/software-compiled-on-EPYC-7313/ASTER/bin/astral-hybrid"

# Define input files
DATADIR="$PROJDIR/data/mirarab2015astral2-extsim/$MODL/$REPL"
STRE_TRUE="s_tree.trees"
GTRE="estimatedgenetre"
ROPTS="-n 0.0 -x 1.0 -d 0.0"
if [ $SUPP == "abayes" ]; then
    ROPTS="-n 0.333 -x 1.0 -d 0.0"
    GTRE="$GTRE.abayes"
fi

# Do work
cd $DATADIR

GTRE_FILE="${GTRE}.${NGEN}"
if [ ! -e $GTRE_FILE ]; then
    head -n${NGEN} $GTRE > $GTRE_FILE
fi

MYMTHDS=( "aster_h_t16" \
	  "wtreeqmc_wh_n2" )

MYSTRE="true_stree_${SUPP}_${NGEN}gen"
if [ ! -e $MYSTRE.tre ]; then
    if [ ! -e $STRE_TRUE ]; then
        echo "$DATADIR/$STRE_TRUE is missing!"
    else
        if [ -z $(grep ";" $STRE_TRUE) ]; then
            echo "$DATADIR/$STRE_TRUE is empty!"
        else
            # Run weighted ASTRAL
            $ASTERH $ROPTS \
                    -u 2 -t 16 \
                    -i $GTRE_FILE \
                    -c $STRE_TRUE \
                    -o $MYSTRE-ah-annotated.tre \
                    &> $MYSTRE-ah-annotated.log

            # Write to CSV file
            MYQSCR="$(grep "Score:" $MYSTRE-ah-annotated.log | awk '{print $2}')"
            echo "$MYMODL,TRUE,$MYQSCR" > ${MYSTRE}_ah_quartet_score.csv
        fi
    fi
fi

for MYMTHD in ${MYMTHDS[@]}; do
    MYSTRE="${MYMTHD}_${SUPP}_${NGEN}gen"
    if [ ! -e $MYSTRE-ah-annotated.tre ]; then
	if [ ! -e $MYSTRE.tre ]; then
            echo "$DATADIR/$MYSTRE.tre is missing!"
        else
            if [ -z $(grep ";" $MYSTRE.tre) ]; then
                echo "$DATADIR/$MYSTRE.tre is empty!"
            else
                # Run weighted ASTRAL
                $ASTERH $ROPTS \
                    -u 2 -t 16 \
                    -i $GTRE_FILE \
                    -c $MYSTRE.tre \
                    -o $MYSTRE-ah-annotated.tre \
                    &> $MYSTRE-ah-annotated.log

                # Write to CSV file
                MYQSCR="$(grep "Score:" $MYSTRE-ah-annotated.log | awk '{print $2}')"
                echo "$MYMODL,$MYMTHD,$MYQSCR" > ${MYSTRE}_ah_quartet_score.csv
            fi
        fi
    fi
done

