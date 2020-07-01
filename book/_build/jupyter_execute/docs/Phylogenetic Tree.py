# Phylogenetic Tree Analysis

## What is phylogenetic tree?
Phylogenetics is the study of the evolutionary history and relationships among individuals or groups of organisms

## Key points
- A phylogenetic tree is a diagram that represents evolutionary relationships among organisms. Phylogenetic trees are hypotheses, not definitive facts
- The pattern of branching in a phylogenetic tree reflects how species or other groups evolved from a series   of common ancestors.
- In trees, two species are more related if they have a more recent common ancestor and less related if they have a less recent common ancestor.
- Phylogenetic trees can be drawn in various equivalent styles. Rotating a tree about its branch points doesn't change the information it carries.

## Types of phylogenetic tree
- Rooted : Rooted tree directed to a unique node.
- Unrooted : Unrooted tree shows the relatedness of the leaves
  without assuming ancestry at all 

## Anatomy of a tree
- Root: origin of evolution
- Leaves: current organisms, species, or genomic
  sequence 
- Branches: relationship between organisms, species,
    or genomic sequence 
- Branch length: evolutionary time 
![Phylo](../figs/phylo1.png)

## Most recent common ancestor
![Phylo](../figs/phylo2.png)

## Which species are more related?
![Phylo](../figs/phylo3.png)
## How to read phylogenetic tree?
### Analysis-1
![Phylo](../figs/phylo4.png)
### Analysis-1
![Phylo](../figs/phylo5.png)
## Steps of constructing phylogenetic tree
- Step 1: Acquiring the Sequences
    - DNA Sequences
    - Protein Sequences
- Step 2: Multiple sequence alignment
    - MEGA
    - ClustalW
- Step 3: Model/Algorithms selection
    - UPGAMA
    - NJ
- Step 4: Phylogenetic tree construction
    - Distance based methods
    - Probabilistic methods
    - Maximum Parsimony (MP) methods
- Step 5: Evaluation/Analysis of tree
    - Statistical analysis

# Implementation in Python

## DistanceTreeConstructor
- The DistanceTreeConstructor has two algorithms:
    - UPGMA (Unweighted Pair Group Method with Arithmetic Mean)
    - NJ (Neighbor Joining)
- Both algorithms construct trees based on a distance matrix. 
- So before using these algorithms, you have to calculate distance matrix from a multiple sequence alignment object by using **DistanceCalculator**

# Essential imports 
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import AlignIO

# Read sequence(MSA)
seq = AlignIO.read("../data/msa.phy", "phylip")
print(seq)

## Distance Matrix Calculation

# Calculate distance matrix 
cal = DistanceCalculator('identity')
dm = cal.get_distance(seq)
print(dm)

## Tree Construction: UPGMA 

tree_constructor = DistanceTreeConstructor()
tree = tree_constructor.upgma(dm)
Phylo.draw(tree)
Phylo.draw_ascii(tree)

## Tree Construction: NJ

tree_constructor = DistanceTreeConstructor()
tree = tree_constructor.nj(dm)
Phylo.draw(tree)
Phylo.draw_ascii(tree)

## References 
- https://www.khanacademy.org/science/high-school-biology/hs-evolution/hs-phylogeny/a/phylogenetic-trees
- https://biopython.org/wiki/Phylo
- https://www.pellegrini.mcdb.ucla.edu/wp-content/uploads/sites/21/2017/07/week-3c-Phylogenetic_Tree_ConstructionMai-copy.pdf
- https://academic.oup.com/mbe/article/30/5/1229/992850