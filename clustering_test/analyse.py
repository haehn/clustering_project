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
parser.add_argument('-c', '--method', type=str)
parser.add_argument('-i', '--input', type=fpath)
parser.add_argument('-p', '--program', help='which program to use. can be used in conjunction with (-m, --model), (-n, --ncat) and (-data, --datatype)', type=str, default='bionj')
parser.add_argument('-m', '--model', help='which model to use in phylogenetic inference', default=None)
parser.add_argument('-n', '--ncat', help='number of categories for gamma distributed rates', default=1)
parser.add_argument('-data', '--datatype', help='datatype', default=None)

args = vars(parser.parse_args())

dist = args['distance']
method = args['method']
path = args['input']
program = args['program']
model = args['model']
ncat = args['ncat']
datatype = args['datatype']
tmp = os.environ['TEMPORARY_DIRECTORY']
gtp_path = os.environ['GTP_PATH']


seq = cPickle.load(open('{0}/seq.pickle'.format(path)))
seq.put_partitions(dist, method, 4, tmpdir=tmp, gtp_path=gtp_path)
seq.put_clusters()
seq.put_cluster_trees(program=program, model=model, datatype=datatype, ncat=ncat, overwrite=False)
result = seq.get_clusters()[(dist,method,4)]
partition = seq.get_partitions()[(dist,method,4)]
truth = seq.get_partitions()['true']
varinf = seq.clustering.variation_of_information(partition, truth)
with open('{0}/{1}{2}.txt'.format(path,dist,method),'w') as file:
	file.write('''Distance:\t{0}
Method:\t{1}
Score:\t{2}
Varinf:\t{3}
'''.format(dist,method,result.score,varinf))
