#!/bin/bash

exit

GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
OUTDIR="$PROJDIR/data/wu2024genomes/species-trees"
WTREEQMC="$PROJDIR/software/TREE-QMC/tree-qmc"

cd $OUTDIR

#diff treeA_wastrid_support_none_nonefilter.tre wastrid_support_cds_and_introns_raxml_abayes.tre 
#diff treeB_wtreeqmc_hybrid_length.tre wtreeqmc_hybrid_n2_cds_and_introns_raxml_abayes.tre
#diff treeC_wtreeqmc_support_none_nonefilter.tre wtreeqmc_support_n2_cds_and_introns_raxml_abayes.tre
#diff treeD_waster_hybrid.tre waster_hybrid_t16_cds_and_introns_raxml_abayes.tre
#diff treeE_waster_length.tre waster_length_t16_cds_and_introns_raxml_abayes.tre 
#diff treeF_waster_support_none.tre waster_support_t16_cds_and_introns_raxml_abayes.tre
#treeG_waster_nonefilter.tre waster_none_t16_cds_and_introns_raxml_abayes_treeshrink.tre

TREES=( "treeA_wastrid_support_none_nonefilter.tre" \
        "treeB_wtreeqmc_hybrid_length.tre" \
        "treeC_wtreeqmc_support_none_nonefilter.tre" \
        "treeD_waster_hybrid.tre" \
        "treeE_waster_length.tre" \
        "treeF_waster_support_none.tre" \
        "treeG_waster_nonefilter.tre" )

for STRE in ${TREES[@]}; do

$WTREEQMC \
    --rootonly $STRE \
    --root allgtrmss \
    -i ../cds_and_introns_raxml_abayes.tre \
    -o rooted_${STRE}

python3 ../compare_two_trees.py -t1 $STRE -t2 rooted_${STRE}

done

