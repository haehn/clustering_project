#!/usr/bin/python
# -*- coding: utf-8 -*-

from sequence_collection import SequenceCollection
import cPickle
import time
import matplotlib.pyplot as plt

indir = '/Users/kgori/git/kevin/yeast_data/MSA'

def add_to_plot(SecCol, metric, link, *args):
    cl = SecCol.get_clusters()
    x,y = [],[]
    for k in sorted(cl):
        if metric in k and link in k:
            y.append(cl[k].score)
            x.append(k[-1])
    plt.plot(x,y, *args)

print 'test directory = ', indir

col = cPickle.load(file('yeast.pickle'))

randomised_collections = []

for i in range(100):
    rand = SequenceCollection(records=col.get_randomised_alignments(),datatype='dna')
    rand.put_trees_parallel()
    rand.put_partitions(metrics=['sym'],linkages=['ward'],nclasses=[1,2,3,4,5,6,7,8,9,10,11,12])
    rand.put_clusters()
    rand.put_cluster_trees()

cPickle.dump(randomised_collections, file('randlist.pickle','w'))   

lines = ['r--','y--','g--','c--','m--','k--']
fig = plt.figure()
add_to_plot(col,'sym','ward')
for i in range(len(randomised_collections)):
    add_to_plot(randomised_collections[i], 'sym', 'ward', lines[i%len(lines)])
fig.show()

