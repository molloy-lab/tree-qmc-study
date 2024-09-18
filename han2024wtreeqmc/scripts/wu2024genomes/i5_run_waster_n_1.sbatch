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

# Define directories
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="/fs/cbcb-lab/ekmolloy/ekmolloy/tree-qmc-study/han2024wtreeqmc"
GENEDIR="$PROJDIR/data/wu2024genomes/gene-trees-abayes"
OUTDIR="$PROJDIR/data/wu2024genomes/species-trees"
ASTERH="$GROUPDIR/group/software-compiled-on-EPYC-7313/ASTER/bin/astral-hybrid"
# --mode	Integer	1	1: hybrid weighting, 2: support only, 3: length only, 4: unweighted

THREAD=1

cd $OUTDIR

# ASTER w/ no weighting scheme
MTHD="waster_none_t${THREAD}"
GTRE="cds_and_introns_raxml_abayes"
STRE="${MTHD}_${GTRE}"
if [ ! -e $STRE.tre ]; then
    TIME="$(time ($ASTERH -t $THREAD -u 0 --mode 4 -i $GENEDIR/${GTRE}.tre -o $STRE.tre &> $STRE.log) 2>&1 1>/dev/null)"

    uname -a > ${STRE}_node_info.csv
    NODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

    SECS="$(echo $TIME | awk '{print $2","$4","$6}')"
    echo "$GTRE,$MTHD,$NODE,$SECS" > ${STRE}_runtime.csv
fi

