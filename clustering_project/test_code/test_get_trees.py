#!/usr/bin/python
# -*- coding: utf-8 -*-


from sequence_collection import SequenceCollection
import time

indir = '/Users/kgori/git/kevin/data/real_data/yeast_data/MSA'

print 'test directory = ', indir

load_start = time.time()
print 'loading sequences (parallel)'
col = SequenceCollection(indir, datatype='dna')
print col
load_end = time.time()

col.put_trees_parallel()
col.put_partitions(metrics=['sym','geodesic'], linkages=['ward'], nclasses=[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16])
col.put_clusters()
col.put_cluster_trees_parallel()

timings = [
load_end - load_start
]

print 'time = {0:.3f}'.format(*timings)
