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

args = vars(parser.parse_args())
outdir = args['directory']

seq = SequenceCollection('{0}/dna_alignments'.format(outdir), datatype='dna', helper=os.environ['DARWINHELPER'])
seq.put_trees(program='bionj')
print 'doing geodesic distance matrices'
seq.put_distance_matrices('geodesic')
print 'doing euc distance matrices'
seq.put_distance_matrices('euc')
print 'doing sym distance matrices'
seq.put_distance_matrices('sym')

print 'pickling into {0}'.format(outdir)
cPickle.dump(seq, file('{0}/seq.pickle'.format(outdir),'w'))
