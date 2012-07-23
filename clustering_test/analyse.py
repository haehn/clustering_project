#!/usr/bin/env python

import os
import cPickle
import argparse

def fpath(s):
    """
    Helper function used when passing filepath arguments with argparse module.
        Trims all '/' characters from the end of the path string.
    """
    while s.endswith('/'):
        s = s[:-1]
    return s

parser = argparse.ArgumentParser(prog='analyse.py')
parser.add_argument('-d', '--distance', type=str)
parser.add_argument('-m', '--method', type=str)
parser.add_argument('-i', '--input', type=fpath)

args = vars(parser.parse_args())

dist = args['distance']
method = args['method']
pickle = args['input']

seq = cPickle.load(file('{0}/seq.pickle'.format(pickle)))
seq.put_partitions(dist, method, 4)
seq.put_clusters()
seq.put_cluster_trees(program='bionj')
result = seq.get_clusters()[(dist,method,4)]
partition = seq.get_partitions()[(dist,method,4)]
with open('{0}/{1}{2}.txt'.format(pickle,dist,method),'w') as file:
	file.write('{0}\n{1}\n'.format)
