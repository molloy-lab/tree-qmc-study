#!/bin/bash
#SBATCH --qos=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64gb
#SBATCH --account=cbcb
#SBATCH --partition=cbcb
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --constraint=EPYC-7313

GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="/fs/cbcb-lab/ekmolloy/ekmolloy/tree-qmc-study/han2024wtreeqmc"
GENEDIR="$PROJDIR/data/wu2024genomes/gene-trees-abayes"
OUTDIR="$PROJDIR/data/wu2024genomes/species-trees"
ASTEROID="$GROUPDIR/group/software-compiled-on-EPYC-7313/Asteroid-nompi/Asteroid/build/bin/asteroid"

cd $OUTDIR

MTHD="asteroid"

# Asteroid
GTRE="cds_and_introns_raxml_abayes"
STRE="${MTHD}_${GTRE}"
if [ ! -e $STRE ]; then
    TIME="$(time ($ASTEROID -i $GENEDIR/${GTRE}.tre -p $STRE.tre &> $STRE.log) 2>&1 1>/dev/null)"

    uname -a > ${STRE}_node_info.csv
    NODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

    SECS="$(echo $TIME | awk '{print $2","$4","$6}')"
    echo "$GTRE,$MTHD,$NODE,$SECS" > ${STRE}_runtime.csv
fi

# Asteroid after treeshrink taxon filtering
GTRE="cds_and_introns_raxml_abayes_treeshrink"
STRE="${MTHD}_${GTRE}"
echo $STRE.tre
ls ../$GTRE.tre
if [ ! -e $STRE ]; then
    TIME="$(time ($ASTEROID -i $GENEDIR/${GTRE}.tre -p $STRE.tre &> $STRE.log) 2>&1 1>/dev/null)"

    uname -a > ${STRE}_node_info.csv
    NODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

    SECS="$(echo $TIME | awk '{print $2","$4","$6}')"
    echo "$GTRE,$MTHD,$NODE,$SECS" > ${STRE}_runtime.csv
fi

