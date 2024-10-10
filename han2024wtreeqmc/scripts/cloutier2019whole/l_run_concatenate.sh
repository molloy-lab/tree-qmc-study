#!/bin/bash
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
INDIR="$GROUPDIR/group/data/cloutier2019whole/alignments"
OUTDIR="$PROJDIR/data/cloutier2019whole/alignments"

CONCAT="$PROJDIR/tools/seqtools.py -c -f fasta"
KEEP="$PROJDIR/tools/seqtools.py -f fasta -k anoDid,aptHaa,aptMan,aptOwe,aptRow,casCas,cryCin,droNov,eudEle,notPer,rheAme,rhePen,strCam,tinGut"

cd $OUTDIR

# Concatenate introns
OUTPUT="concat_introns.fasta"
if [ ! -e $OUTPUT ]; then
    INPUT=$(ls $INDIR/introns/*.fasta)
    python3 $CONCAT -i $INPUT -o tmp-$OUTPUT
    sed 's/N/-/g' tmp-$OUTPUT | sed 's/dro-ov/droNov/g' > $OUTPUT
fi
if [ ! -e nogal-$OUTPUT ]; then
    python3 $KEEP -i $OUTPUT -o nogal-$OUTPUT
fi


# Concatenate UCEs
OUTPUT="concat_UCEs_plus_bad105.fasta"
if [ ! -e $OUTPUT ]; then
    INPUT=$(ls $INDIR/UCEs/*.fasta)
    python3 $CONCAT -i $INPUT -o tmp-$OUTPUT
    sed 's/N/-/g' tmp-$OUTPUT | sed 's/dro-ov/droNov/g' > $OUTPUT
fi
if [ ! -e nogal-$OUTPUT ]; then
    python3 $KEEP -i $OUTPUT -o nogal-$OUTPUT
fi


# Concatenate good UCEs (no 105 bad UCEs)
OUTPUT="concat_UCEs.fasta"
if [ ! -e $OUTPUT ]; then
    BADLIST=( $(cat $PROJDIR/scripts/cloutier2019whole/list_of_105_UCEs_to_exclude.txt) )
    mkdir -p $INDIR/UCEs/exclude
    for UCE in ${BADLIST[@]}; do
        mv $INDIR/UCEs/$UCE.fasta $INDIR/UCEs/exclude
    done
    INPUT=$(ls $INDIR/UCEs/*.fasta)
    python3 $CONCAT -i $INPUT -o tmp-$OUTPUT
    mv $INDIR/UCEs/exclude/* $INDIR/UCEs
    sed 's/N/-/g' tmp-$OUTPUT | sed 's/dro-ov/droNov/g' > $OUTPUT
fi
if [ ! -e nogal-$OUTPUT ]; then
    python3 $KEEP -i $OUTPUT -o nogal-$OUTPUT
fi


# Concatenate bad UCEs but minus impacted taxa (galGal and tinGut)
OUTPUT="bad105_without_galGal_tinGut.fasta"
if [ ! -e $OUTPUT ]; then
    BADLIST=( $(cat $PROJDIR/scripts/cloutier2019whole/list_of_105_UCEs_to_exclude.txt) )
    mkdir -p $INDIR/UCEs/exclude
    for UCE in ${BADLIST[@]}; do
        mv $INDIR/UCEs/$UCE.fasta $INDIR/UCEs/exclude
    done
    INPUT=$(ls $INDIR/UCEs/exclude/*.fasta)
    python3 $CONCAT -i $INPUT -o tmp-$OUTPUT -k "anoDid,aptHaa,aptMan,aptOwe,aptRow,casCas,cryCin,droNov,eudEle,notPer,rheAme,rhePen,strCam"
    mv $INDIR/UCEs/exclude/* $INDIR/UCEs
    sed 's/N/-/g' tmp-$OUTPUT | sed 's/dro-ov/droNov/g' > $OUTPUT
fi
if [ ! -e nogal-$OUTPUT ]; then
    python3 $KEEP -i $OUTPUT -o nogal-$OUTPUT
fi


# Concatenate good and bad minus impacted taxa (galGal and tinGut)
OUTPUT="concat_UCEs_plus_bad105_without_galGal_tinGut.fasta"
if [ ! -e $OUTPUT ]; then
    python3 $CONCAT -i concat_UCEs.fasta bad105_without_galGal_tinGut.fasta -o tmp-$OUTPUT
    sed 's/N/-/g' tmp-$OUTPUT | sed 's/dro-ov/droNov/g' > $OUTPUT
fi
if [ ! -e nogal-$OUTPUT ]; then
    python3 $KEEP -i $OUTPUT -o nogal-$OUTPUT
fi

