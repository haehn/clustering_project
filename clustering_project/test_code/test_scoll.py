#!/usr/bin/python
# -*- coding: utf-8 -*-

from sequence_record import TCSeqRec
from sequence_collection import SequenceCollection
from tree import Tree
from clustering import Clustering
import cPickle
import time
import os
import copy
import numpy as np
np.set_printoptions(precision=2, linewidth=200)

indir = '/Users/kgori/git/kevin/data/simulated_data/small/MSA'

print 'test directory = ', indir

load_start = time.time()
print 'loading sequences (parallel)'
col = SequenceCollection(indir, datatype='protein')
print col
load_end = time.time()

tcseq_start = time.time()
print 'getting TC trees (sequential)'
col.get_trees(program='treecollection', tmpdir='/tmp')
for rec in col.records:
    print rec.name
    print rec.tree
tcseq_end = time.time()

tcpar_start = time.time()
print 'getting TC trees (parallel)'
col.get_trees_parallel(program='treecollection', tmpdir='/tmp')
for rec in col.records:
    print rec.name
    print rec.tree
tcpar_end = time.time()

raxseq_start = time.time()
print 'getting raxml trees (sequential)'
col.get_trees(program='raxml', tmpdir='/tmp')
for rec in col.records:
    print rec.name
    print rec.tree
raxseq_end = time.time()

raxpar_start = time.time()
print 'getting raxml trees (parallel)'
col.get_trees_parallel(program='raxml', tmpdir='/tmp')
for rec in col.records:
    print rec.name
    print rec.tree
raxpar_end = time.time()

physeq_start = time.time()
print 'getting phyml trees (sequential)'
col.get_trees(program='phyml', tmpdir='/tmp')
for rec in col.records:
    print rec.name
    print rec.tree
physeq_end = time.time()

phypar_start = time.time()
print 'getting phyml trees (parallel)'
col.get_trees_parallel(program='phyml', tmpdir='/tmp')
for rec in col.records:
    print rec.name
    print rec.tree
phypar_end = time.time()

timings = [
    load_end - load_start,
    tcseq_end - tcseq_start,
    tcpar_end - tcpar_start,
    raxseq_end - raxseq_start,
    raxpar_end - raxpar_start,
    physeq_end - physeq_start,
    phypar_end - phypar_start,
    phypar_end - load_start
    ]

print '''
Timings:
loading          = {0:.3f}
Sequential TC    = {1:.3f}
Parallel TC      = {2:.3f}
Sequential raxml = {3:.3f}
Parallel raxml   = {4:.3f}
Sequential phyml = {5:.3f}
Parallel phyml   = {6:.3f}
Total            = {7:.3f}
'''.format(*timings)

# x = SequenceCollection(indir, datatype='dna')
# cPickle.dump(x, file('yeast_x.pickle', 'w'))
# y = copy.deepcopy(x)
# z = copy.deepcopy(x)
# x = cPickle.load(file('yeast_x.pickle'))
# x.get_TC_trees_parallel()

## y.get_raxml_trees_parallel()
## z.get_phyml_trees_parallel()

# for collection in [x]:
#     for metric in ['sym','rf','euc','geodesic']:
#         collection.get_distance_matrix(metric, normalise=False)
#         for linkage in ['single','complete','ward']:
#             for n in range(2,9):
#                 collection.get_partitions(metric, linkage, n)

# cPickle.dump(x, file('yeast_x.pickle', 'w'))
# print 'Getting clusters'
# x.get_clusters()
# print 'Getting cluster trees in parallel'
# x.get_cluster_trees_parallel()
# print 'checking everything is where it should be...'
# for k in x.clustering.clusters:
#     for cl in x.clustering.clusters[k]:
#         print cl['cluster'].tree

## e = SequenceCollection(indir, datatype='protein')
# t_pre_parallel = time.time()
# y.get_TC_trees_parallel()

# t_post_parallel = time.time()
# y.get_TC_trees()
# t_post_sequential = time.time()
# print 'timings: parallel: {0:.3f}s; sequential: {1:.3f}s'.format(t_post_parallel - t_pre_parallel, t_post_sequential - t_pre_parallel)
# y.get_TC_trees()
# t1 = time.time()
# y.get_distance_matrix('sym', normalise=False)
# t2 = time.time()
# y.get_distance_matrix('rf', normalise=False)
# t3 = time.time()
# y.get_distance_matrix('euc', normalise=False)
# t4 = time.time()
# y.get_distance_matrix('geodesic')
# t5 = time.time()
# for k in y.clustering.distance_matrices:
#     print k
#     print y.clustering.distance_matrices[k]
#     for l in ['single', 'complete', 'ward']:
#         y.get_clustering(k, l, 4)
# t6 = time.time()

# print 'timings: sym: {0}s; rf: {1}s; euc:{2}s; geo:{3}s; clustering:{4}s; total: {5}s'.format(
#     t2 - t1,
#     t3 - t2,
#     t4 - t3,
#     t5 - t4,
#     t6 - t5,
#     t6 - t1,
#     )

# d.clustering.concatenate_records(('geodesic','ward',4),d.records)

# for k in d.clustering.distance_matrices:
#     for l in ['single', 'complete', 'ward']:
#         print k, l
#         p = d.plot_clustering(k, l, 4, show=True)
#         p.show()
# d=cPickle.load(file('d_para.pickle'))
# print d.clustering

# cPickle.dump(x, file('yeast_tc.pickle', 'w'))
# cPickle.dump(y, file('yeast_rax.pickle', 'w'))
# cPickle.dump(z, file('yeast_phy.pickle', 'w'))
