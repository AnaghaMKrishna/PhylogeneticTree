# PhylogeneticTree
Create a phylogenetic tree that represents evolutionary relationships among organisms. The pattern of branching in a phylogenetic tree reflects how species or other groups evolved from a series of common ancestors.

**phylogenetic_tree.py** generates phylogenetic tree in newick format by using UPGMA(Unweighted Pair Group Method with Arithmetic mean) to identify hierarchical relationships among given organisms.

**tree_multi_primer.py** extends the above implementation by considering multiple loci, thus providing higher resolution, for evaluating relatedness among organisms. 

**phyl_tree_with_tools/phyl_tree_tool.py** creates phylogenetic tree in newick format considering both single locus and multiple loci using Muscle for Multi-Sequence Alignment(MSA) and FastTree for generating maximum-likelihood tree.