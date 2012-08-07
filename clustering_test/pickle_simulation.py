#!/usr/bin/env python

import os
import cPickle
import argparse
from sequence_collection import SequenceCollection

def fpath(s):
    """
    Helper function used when passing filepath arguments with argparse module.
        Trims all '/' characters from the end of the path string.
    """
    while s.endswith('/'):
        s = s[:-1]
    return s

parser = argparse.ArgumentParser(prog='pickle.py')
parser.add_argument('-d', '--directory', help='output directory', type=fpath, default='.')
parser.add_argument('-g', '--geodesic', help='path to gtp.jar', type=fpath, default='.')
parser.add_argument('-p', '--program', help='which program to use. can be used in conjunction with (-m, --model), (-n, --ncat) and (-data, --datatype)', type=str, default='bionj')
parser.add_argument('-m', '--model', help='which model to use in phylogenetic inference', default=None)
parser.add_argument('-n', '--ncat', help='number of categories for gamma distributed rates', default=1)
parser.add_argument('-t', '--tmpdir', help='temporary directory', type=fpath, default='/tmp')
parser.add_argument('-data', '--datatype', help='datatype', default=None)

args = vars(parser.parse_args())
outdir = args['directory']
program = args['program']
model = args['model']
ncat = args['ncat']
datatype = args['datatype']
gtp_path = os.environ['GTP_PATH']
tmpdir = os.environ['TEMPORARY_DIRECTORY']
tmpdir = args['tmpdir']
print 'Reading alignments into SequenceRecord object'
seq = SequenceCollection('{0}/dna_alignments'.format(outdir), datatype='dna', tmpdir=tmpdir, helper=os.environ['DARWINHELPER'])
print 'Calculating trees'
print program,model,datatype,ncat,tmpdir
seq.put_trees_parallel(program=program, model=model, datatype=datatype, ncat=ncat, tmpdir=tmpdir)
print 'doing geodesic distance matrices'
seq.put_distance_matrices('geo', gtp_path=gtp_path, tmpdir=tmpdir)
print 'doing euc distance matrices'
seq.put_distance_matrices('euc')
print 'doing sym distance matrices'
seq.put_distance_matrices('sym')

print 'getting score for true clustering'
with open('{0}/treedistances.txt'.format(outdir)) as file:
	T = file.readline().rstrip().split('\t')[1][1:-1].split(', ')
seq.clustering.partitions['true'] = T
seq.put_clusters()
seq.put_cluster_trees_parallel(program=program, model=model, datatype=datatype, ncat=ncat, tmpdir=tmpdir)
r = seq.get_clusters()['true']
print 'score: ', r.score
with open('{0}/treedistances.txt'.format(outdir),'a') as file:
    file.write('True score: {0}\n'.format(r.score))

print 'pickling into {0}'.format(outdir)
cPickle.dump(seq, open('{0}/seq.pickle'.format(outdir),'w'))
