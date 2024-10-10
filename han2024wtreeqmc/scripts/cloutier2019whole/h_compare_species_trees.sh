#!/bin/bash

GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
DATADIR="$PROJDIR/data/cloutier2019whole/species-trees-analysis"
SIMMONSDIR="$GROUPDIR/group/data/simmons2022gene/Simmons_et_al_trees_and_matrices"
COMPARE="$PROJDIR/tools/compare_two_trees.py"

cd $DATADIR

TA="treeA.tre"
TB="treeB.tre"
TC="treeC.tre"

python3 $COMPARE -t1 $TA -t2 ../species_trees/ExaML/intron_examl_species_tree.tre
python3 $COMPARE -t1 $TC -t2 ../species_trees/ExaML/UCE_examl_species_tree.tre
python3 $COMPARE -t1 $TA -t2 ../species_trees/ExaML/TENT_examl_species_tree.tre

python3 $COMPARE -t1 $TA -t2 $SIMMONSDIR/concatenation/likelihood/plus_Gallus/RAxML_bipartitions.intron_mafft_BS.tre 
python3 $COMPARE -t1 $TC -t2 $SIMMONSDIR/concatenation/likelihood/plus_Gallus/RAxML_bipartitions.UCE_mafft_BS.tre 
python3 $COMPARE -t1 $TA -t2 $SIMMONSDIR/concatenation/likelihood/plus_Gallus/RAxML_bipartitions.TENT_mafft_BS.tre 
python3 $COMPARE -t1 $TC -t2 $SIMMONSDIR/UCE_MAFFT_alignment_homology_errors/concatenation/RAxML_bipartitions.UCE_minus_105_plus_Gallus.tre

# Check the ASTEROID trees
python3 $COMPARE -t1 $TB -t2 asteroid_introns_abayes_gene_trees_sorted.tre.bestTree.newick
python3 $COMPARE -t1 $TB -t2 asteroid_UCEs_minus_bad105_sorted_plus_bad105_sorted_abayes_gene_trees.tre.bestTree.newick
python3 $COMPARE -t1 $TB -t2 asteroid_UCEs_minus_bad105_sorted_plus_bad105_sorted_without_galGal_and_tinGut_abayes_gene_trees.tre.bestTree.newick
python3 $COMPARE -t1 $TB -t2 asteroid_UCEs_minus_bad105_abayes_gene_trees_sorted.tre.bestTree.newick

# Check unweighted ASTRID trees
python3 $COMPARE -t1 $TB -t2 wastrid_none_introns_abayes_gene_trees_sorted.tre
python3 $COMPARE -t1 $TB -t2 wastrid_none_UCEs_minus_bad105_sorted_plus_bad105_sorted_abayes_gene_trees.tre
python3 $COMPARE -t1 $TB -t2 wastrid_none_UCEs_minus_bad105_sorted_plus_bad105_sorted_without_galGal_and_tinGut_abayes_gene_trees.tre
python3 $COMPARE -t1 $TB -t2 wastrid_none_UCEs_minus_bad105_abayes_gene_trees_sorted.tre

# Check support weighted ASTRID tree
python3 $COMPARE -t1 $TB -t2 wastrid_support_introns_abayes_gene_trees_sorted.tre
python3 $COMPARE -t1 $TB -t2 wastrid_support_UCEs_minus_bad105_sorted_plus_bad105_sorted_abayes_gene_trees.tre
python3 $COMPARE -t1 $TB -t2 wastrid_support_UCEs_minus_bad105_sorted_plus_bad105_sorted_without_galGal_and_tinGut_abayes_gene_trees.tre
python3 $COMPARE -t1 $TB -t2 wastrid_support_UCEs_minus_bad105_abayes_gene_trees_sorted.tre

# Check length weighted ASTRID trees
python3 $COMPARE -t1 $TB -t2 wastrid_length_introns_abayes_gene_trees_sorted.tre
python3 $COMPARE -t1 $TB -t2 wastrid_length_UCEs_minus_bad105_sorted_plus_bad105_sorted_abayes_gene_trees.tre
python3 $COMPARE -t1 $TB -t2 wastrid_length_UCEs_minus_bad105_sorted_plus_bad105_sorted_without_galGal_and_tinGut_abayes_gene_trees.tre
python3 $COMPARE -t1 $TB -t2 wastrid_length_UCEs_minus_bad105_abayes_gene_trees_sorted.tre

# Check unweighted TREE-QMC 
python3 $COMPARE -t1 $TB -t2 wtreeqmc_none_n2_introns_abayes_gene_trees_sorted.tre
python3 $COMPARE -t1 $TB -t2 wtreeqmc_none_n2_UCEs_minus_bad105_sorted_plus_bad105_sorted_abayes_gene_trees.tre
python3 $COMPARE -t1 $TA -t2 wtreeqmc_none_n2_UCEs_minus_bad105_sorted_plus_bad105_sorted_without_galGal_and_tinGut_abayes_gene_trees.tre
python3 $COMPARE -t1 $TA -t2 wtreeqmc_none_n2_UCEs_minus_bad105_abayes_gene_trees_sorted.tre

# Check support-weighted TREE-QMC 
python3 $COMPARE -t1 $TB -t2 wtreeqmc_support_n2_introns_abayes_gene_trees_sorted.tre
python3 $COMPARE -t1 $TB -t2 wtreeqmc_support_n2_UCEs_minus_bad105_sorted_plus_bad105_sorted_abayes_gene_trees.tre
python3 $COMPARE -t1 $TA -t2 wtreeqmc_support_n2_UCEs_minus_bad105_sorted_plus_bad105_sorted_without_galGal_and_tinGut_abayes_gene_trees.tre
python3 $COMPARE -t1 $TA -t2 wtreeqmc_support_n2_UCEs_minus_bad105_abayes_gene_trees_sorted.tre

# Check length-weighted TREE-QMC 
python3 $COMPARE -t1 $TA -t2 wtreeqmc_length_n2_introns_abayes_gene_trees_sorted.tre
python3 $COMPARE -t1 $TA -t2 wtreeqmc_length_n2_UCEs_minus_bad105_sorted_plus_bad105_sorted_abayes_gene_trees.tre
python3 $COMPARE -t1 $TA -t2 wtreeqmc_length_n2_UCEs_minus_bad105_sorted_plus_bad105_sorted_without_galGal_and_tinGut_abayes_gene_trees.tre
python3 $COMPARE -t1 $TA -t2 wtreeqmc_length_n2_UCEs_minus_bad105_abayes_gene_trees_sorted.tre

# Check hybrid-weighted TREE-QMC 
python3 $COMPARE -t1 $TA -t2 wtreeqmc_hybrid_n2_introns_abayes_gene_trees_sorted.tre
python3 $COMPARE -t1 $TA -t2 wtreeqmc_hybrid_n2_UCEs_minus_bad105_sorted_plus_bad105_sorted_abayes_gene_trees.tre
python3 $COMPARE -t1 $TA -t2 wtreeqmc_hybrid_n2_UCEs_minus_bad105_sorted_plus_bad105_sorted_without_galGal_and_tinGut_abayes_gene_trees.tre
python3 $COMPARE -t1 $TA -t2 wtreeqmc_hybrid_n2_UCEs_minus_bad105_abayes_gene_trees_sorted.tre

# Check ASTER
python3 $COMPARE -t1 $TB -t2 waster_none_introns_abayes_gene_trees_sorted.tre
python3 $COMPARE -t1 $TB -t2 waster_none_UCEs_minus_bad105_sorted_plus_bad105_sorted_abayes_gene_trees.tre
python3 $COMPARE -t1 $TA -t2 waster_none_UCEs_minus_bad105_sorted_plus_bad105_sorted_without_galGal_and_tinGut_abayes_gene_trees.tre
python3 $COMPARE -t1 $TA -t2 waster_none_UCEs_minus_bad105_abayes_gene_trees_sorted.tre

# Check support-weighted ASTER
python3 $COMPARE -t1 $TB -t2 waster_support_introns_abayes_gene_trees_sorted.tre
python3 $COMPARE -t1 $TB -t2 waster_support_UCEs_minus_bad105_sorted_plus_bad105_sorted_abayes_gene_trees.tre
python3 $COMPARE -t1 $TA -t2 waster_support_UCEs_minus_bad105_sorted_plus_bad105_sorted_without_galGal_and_tinGut_abayes_gene_trees.tre
python3 $COMPARE -t1 $TA -t2 waster_support_UCEs_minus_bad105_abayes_gene_trees_sorted.tre

# Check length-weighted ASTER
python3 $COMPARE -t1 $TA -t2 waster_length_introns_abayes_gene_trees_sorted.tre
python3 $COMPARE -t1 $TA -t2 waster_length_UCEs_minus_bad105_sorted_plus_bad105_sorted_abayes_gene_trees.tre
python3 $COMPARE -t1 $TA -t2 waster_length_UCEs_minus_bad105_sorted_plus_bad105_sorted_without_galGal_and_tinGut_abayes_gene_trees.tre
python3 $COMPARE -t1 $TA -t2 waster_length_UCEs_minus_bad105_abayes_gene_trees_sorted.tre

# Check hybrid-weighted ASTER
python3 $COMPARE -t1 $TA -t2 waster_hybrid_introns_abayes_gene_trees_sorted.tre
python3 $COMPARE -t1 $TA -t2 waster_hybrid_UCEs_minus_bad105_sorted_plus_bad105_sorted_abayes_gene_trees.tre
python3 $COMPARE -t1 $TA -t2 waster_hybrid_UCEs_minus_bad105_sorted_plus_bad105_sorted_without_galGal_and_tinGut_abayes_gene_trees.tre
python3 $COMPARE -t1 $TA -t2 waster_hybrid_UCEs_minus_bad105_abayes_gene_trees_sorted.tre

