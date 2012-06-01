clustering_project
==================

Tidying up code for the ToL clustering project before pushing to ethz repo

The code in this repository is a python module designed to take a collection 
of sequence alignments (of genes), infer their phylogenetic trees, and
perform hierarchical clustering based on a distance matrix of their tree 
topologies.

The purpose of this is to establish whether there is any underlying structure
to the data

Dependencies:

Python:
numpy
scipy
dendropy

External:
Darwin
RAxML
PhyML
TreeCollection
