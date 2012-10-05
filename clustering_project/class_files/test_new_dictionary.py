#!/usr/bin/env python

print 'importing...'
from sequence_collection import SequenceCollection
from distance_matrix import DistanceMatrix
from partition import Partition
from clustering import Clustering
import sys
import numpy as np
import pylab
print 'done.'


def print_dict(d):
    for k in sorted(d):
        print d


np.set_printoptions(linewidth=200, precision=3)
sc = SequenceCollection(
    '/Users/kgori/scratch/chk/aa_alignments/',
    get_distances=True,
    file_format='phylip',
    helper='/Users/kgori/git/kevin/clustering_project/class_files/DV_wrapper.drw'
        ,
    parallel_load=True,
    gtp_path='/Users/kgori/git/kevin/clustering_project/class_files/',
    tmpdir='/tmp',
    datatype='protein',
    )

sc.put_trees_parallel(program='bionj')
sc.put_partitions(['rf', 'wrf', 'geo', 'euc'], [
    'kmedoids',
    'spectral',
    'MDS',
    'single',
    'complete',
    'average',
    'ward',
    ], [3,4,5])

sc.clusters_to_partitions['true'] = (
    1,
    1,
    1,
    1,
    1,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    )
sc.partitions[(
    1,
    1,
    1,
    1,
    1,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    )] = Partition((
    1,
    1,
    1,
    1,
    1,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    ),
    )
sc.concatenate_records()
sc.put_cluster_trees_parallel(program='bionj')

# print sc.get_trees()
# print trees
# print trees == sc.get_trees()

# sc.put_distance_matrices(['rf','wrf','geo','euc'])
# print sc.distance_matrices

# sc.put_partition('rf','MDS',4)
# print sc.clusters_to_partitions
# print sc.partitions
# sys.exit()

# C = Clustering()
# dm = sc.distance_matrices['rf']
# kmed = C.run_clustering(dm,4,'kmedoids')
# spec = C.run_clustering(dm,4,'spectral')
# mdsc = C.run_clustering(dm,4,'MDS')
# sing = C.run_clustering(dm,4,'single')
# comp = C.run_clustering(dm,4,'complete')
# aver = C.run_clustering(dm,4,'average')
# ward = C.run_clustering(dm,4,'ward')
# dbsc = C.run_clustering(dm,4,'dbscan')

# for x in [kmed,spec,mdsc,sing,comp,aver,ward,dbsc]:
#     print x

# print sys.getsizeof(sc)
