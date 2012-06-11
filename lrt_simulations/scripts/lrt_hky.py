#!/usr/bin/python
# -*- coding: utf-8 -*-

from sequence_collection import SequenceCollection
import cPickle
import sys
import os

indir = sys.argv[1]
datatype = sys.argv[2]
helper = os.environ['DARWINHELPER']
tmpdir = os.environ['TEMPORARY_DIRECTORY']
print tmpdir
print 'setting test directory = ', indir

print 'loading sequences...'
col = SequenceCollection(indir, datatype=datatype, helper=helper, tmpdir=tmpdir)

print 'getting trees...'
col.put_trees_parallel(program='phyml', model='HKY85', datatype='nt',
                       ncat=1, tmpdir=tmpdir)

print 'getting partitions...'
col.put_partitions(metrics=['sym'], linkages=['ward'],
                   nclasses=range(1, 21))
col.put_clusters()
col.put_cluster_trees_parallel(program='phyml', model='HKY85',
                               datatype='nt', ncat=1, tmpdir=tmpdir)

cPickle.dump(col, file('yeast_hky.pickle', 'w'))

