#!/usr/bin/python
# -*- coding: utf-8 -*-

from sequence_collection import SequenceCollection
import cPickle
import time

indir = '/Users/kgori/git/kevin/data/simulated_data/small/MSA'

print 'test directory = ', indir

load_start = time.time()
print 'loading sequences (parallel)'
col = SequenceCollection(indir, datatype='protein')
print col
load_end = time.time()

tcpar_start = time.time()
print 'putting TC trees (parallel)'
col.put_trees_parallel(program='treecollection', tmpdir='/tmp')
for rec in col.records:
    print rec.name
    print rec.tree
tcpar_end = time.time()

par_start = time.time()
print 'Putting partitions'
col.put_partitions(metrics=['sym','euc'], linkages=['ward','single'], nclasses = [3,4,5,6])
print col.get_partitions()
par_end = time.time()

col.put_clusters()

tccls_start = time.time()
# print 'putting TC cluster trees (sequential)'
# col.put_cluster_trees(program='treecollection', tmpdir='/tmp')
# for cl in col.get_cluster_records():
#     print cl.name
#     print cl.tree
tccls_end = time.time()

tcclp_start = time.time()
print 'putting TC cluster trees (parallel)'
col.put_cluster_trees_parallel(program='treecollection',
        tmpdir='/tmp')
cl_recs = col.get_cluster_records()
for cl in cl_recs:
    print cl.name
    print cl.tree
tcclp_end = time.time()

timings = [
    load_end - load_start,
    tcpar_end - tcpar_start,
    par_end - par_start,
    tccls_end - tccls_start,
    tcclp_end - tcclp_start,
    tcclp_end - load_start,
    ]

print '''
Timings:
loading          = {0:.3f}
Parallel TC      = {1:.3f}
Linkages         = {2:.3f}
Sequential TCcl  = {3:.3f}
Parallel TCcl    = {4:.3f}
Total            = {5:.3f}
'''.format(*timings)
