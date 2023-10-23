#!/bin/bash

exit

# Define directories
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJECTDIR="$GROUPDIR/ekmolloy/species-to-tumor-study"

# Define input path
INPATH="$GROUPDIR/group/data/kizilkale2022fast"

# Define output path
OUTPATH="$PROJECTDIR/data/kizilkale2022fast"

for FNUM in `seq 1 10`; do
   cd $INPATH/Extended_Data_Figure_${FNUM}/noisy 
   NSIMS=( $(ls *.SC | sed 's/\.SC//g') )
   for NSIM in ${NSIMS[@]}; do
       TMAT="$INPATH/Extended_Data_Figure_${FNUM}/ground/$NSIM.SC.before_FP_FN_NA"
       NMAT="$INPATH/Extended_Data_Figure_${FNUM}/noisy/$NSIM.SC"

       if [ -e $TMAT ] && [ -e $NMAT ]; then
           OUTDIR="$OUTPATH/${NSIM}"
           if [ ! -d $OUTDIR ]; then
               mkdir $OUTDIR
               cp $TMAT $OUTDIR/ground.CFMatrix
               cp $NMAT $OUTDIR/noisy.CFMatrix
           fi
           echo "$FNUM,$NSIM" >> $OUTPATH/map-extended-figure-to-simulation.csv
       else
           echo "Something wrong with $NMAT"
           #cp $INPATH/Extended_Data_Figure_${FNUM}/ground/$NSIM.ground.CFMatrix $TMAT
       fi
   done
done

