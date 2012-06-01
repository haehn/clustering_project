#!/usr/bin/python
# -*- coding: utf-8 -*-

from sequence_collection import SequenceCollection
import cPickle
import sys
import os

indir = sys.argv[1]
datatype = sys.argv[2]
index = sys.argv[3]
helper = os.environ['DARWINHELPER']
tmpdir = os.environ['TEMPORARY_DIRECTORY']
print tmpdir


def add_to_plot(SecCol, metric, link, *args):
    cl = SecCol.get_clusters()
    x,y = [],[]
    for k in sorted(cl):
        if metric in k and link in k:
            y.append(cl[k].score)
            x.append(k[-1])
    return (x,y, args)
print 'setting test directory = ', indir
"""
print 'loading sequences...'
col = SequenceCollection(indir, helper=helper, tmpdir=tmpdir)

print 'getting trees...'
col.put_trees_parallel(program='phyml',tmpdir=tmpdir)

print 'getting partitions...'
col.put_partitions(metrics=['sym'],linkages=['ward'],nclasses=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17])
#print 'randomizing bytes...'
col.put_clusters()

#print 'immanentizing the eschaton...'
col.put_cluster_trees_parallel(program='phyml',tmpdir=tmpdir)
"""
plottable = []
#plottable.append(add_to_plot(col,'sym','ward'))
col = cPickle.load(file('col.pickle'))
#print 'whipping into frenzy...'
for i in range(1):
    r = SequenceCollection(records=col.get_randomised_alignments(), helper=helper, tmpdir=tmpdir)
    r.put_trees_parallel(program='phyml',tmpdir=tmpdir)
    r.put_partitions(metrics=['sym'],linkages=['ward'],nclasses=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17])
    r.put_clusters()
    r.put_cluster_trees_parallel(program='phyml',tmpdir=tmpdir)   
    plottable.append(add_to_plot(r, 'sym', 'ward'))
    
cPickle.dump(plottable, file('plottable{0}.pickle'.format(index),'w'))
#cPickle.dump(col, file('col.pickle','w'))

