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
GENEDIR="$GROUPDIR/group/data/wu2024genomes/1_besttrees_cds5756"
OUTDIR="$PROJDIR/data/wu2024genomes/gene-trees-abayes"
WASTRID="$GROUPDIR/group/software/wastrid"

cd $OUTDIR

# ASTRID w/ support weighting scheme
MTHD="wastrid_support"
GTRE="cds_and_introns_raxml_abayes"
STRE="${MTHD}_${GTRE}"
if [ ! -e $STRE.tre ]; then
    TIME="$(time ($WASTRID -m support -b 0.333-1 -i $GENEDIR/${GTRE}.tre -o $STRE.tre &> $STRE.log) 2>&1 1>/dev/null)"

    uname -a > ${STRE}_node_info.csv
    NODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

    SECS="$(echo $TIME | awk '{print $2","$4","$6}')"
    echo "$GTRE,$MTHD,$NODE,$SECS" > ${STRE}_runtime.csv
fi

# ASTRID w/ length weighting scheme
MTHD="wastrid_length"
GTRE="cds_and_introns_raxml_abayes"
STRE="${MTHD}_${GTRE}"
if [ ! -e $STRE.tre ]; then
    TIME="$(time ($WASTRID -m n-length -i $GENEDIR/${GTRE}.tre -o $STRE.tre &> $STRE.log) 2>&1 1>/dev/null)"

    uname -a > ${STRE}_node_info.csv
    NODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

    SECS="$(echo $TIME | awk '{print $2","$4","$6}')"
    echo "$GTRE,$MTHD,$NODE,$SECS" > ${STRE}_runtime.csv
fi


# ASTRID w/ no weighting scheme
MTHD="wastrid_none"
GTRE="cds_and_introns_raxml_abayes"
STRE="${MTHD}_${GTRE}"
if [ ! -e $STRE ]; then
    TIME="$(time ($WASTRID -m internode -i $GENEDIR/${GTRE}.tre -o $STRE.tre &> $STRE.log) 2>&1 1>/dev/null)"

    uname -a > ${STRE}_node_info.csv
    NODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

    SECS="$(echo $TIME | awk '{print $2","$4","$6}')"
    echo "$GTRE,$MTHD,$NODE,$SECS" > ${STRE}_runtime.csv
fi

# ASTRID w/ no weighting scheme - but treeshrink taxon filtering
MTHD="wastrid_none"
GTRE="cds_and_introns_raxml_abayes_treeshrink"
STRE="${MTHD}_${GTRE}"
echo $STRE.tre
ls ../$GTRE.tre
if [ ! -e $STRE ]; then
    TIME="$(time ($WASTRID -m internode -i $GENEDIR/${GTRE}.tre -o $STRE.tre &> $STRE.log) 2>&1 1>/dev/null)"

    uname -a > ${STRE}_node_info.csv
    NODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

    SECS="$(echo $TIME | awk '{print $2","$4","$6}')"
    echo "$GTRE,$MTHD,$NODE,$SECS" > ${STRE}_runtime.csv
fi

