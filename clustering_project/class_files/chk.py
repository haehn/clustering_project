#!/usr/bin/env python
from distance_matrix import DistanceMatrix
from tree import Tree

trees = [Tree.new_random_coal(10) for _ in range(10)]
print trees
dm = DistanceMatrix(trees)
print dm
dm.get_distance_matrix('wrf')
print dm.matrix
