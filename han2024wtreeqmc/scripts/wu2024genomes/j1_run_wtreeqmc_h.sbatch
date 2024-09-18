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
#SBATCH --time=72:00:00
#SBATCH --constraint=EPYC-7313

lscpu
cat /proc/meminfo

GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
GENEDIR="$PROJDIR/data/wu2024genomes/gene-trees-abayes"
OUTDIR="$PROJDIR/data/wu2024genomes/species-trees"
WTREEQMC="$PROJDIR/software/TREE-QMC/tree-qmc"

cd $OUTDIR

# TREE-QMC w/ hybrid weighting scheme
MTHD="wtreeqmc_hybrid_n2"
GTRE="cds_and_introns_raxml_abayes"
STRE="${MTHD}_${GTRE}"
if [ ! -e $STRE.tre ]; then
    TIME="$(time ($WTREEQMC -w h -n 0.333 -x 1 -d 0.333 -i $GENEDIR/${GTRE}.tre -o $STRE.tre &> $STRE.log) 2>&1 1>/dev/null)"

    uname -a > ${STRE}_node_info.csv
    NODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

    SECS="$(echo $TIME | awk '{print $2","$4","$6}')"
    echo "$GTRE,$MTHD,$NODE,$SECS" > ${STRE}_runtime.csv
fi

