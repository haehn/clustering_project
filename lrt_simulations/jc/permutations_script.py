#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import cPickle
import os
import sys
from sequence_collection import SequenceCollection

# Command-line arguments handled by argparse module
parser = argparse.ArgumentParser(prog='permutations_script.py',
                description='Run random permutations procedure on SequenceCollection'
                )
parser.add_argument(
    '-b',
    '--basefile',
    dest='base',
    help='path to base SequenceCollection object stored as pickle',
    type=str,
    default=''
    )
parser.add_argument(
    '-o',
    '--outfile',
    dest='out',
    help='output file',
    type=str,
    default=''
    )
parser.add_argument(
    '-d',
    '--distance',
    dest='dist',
    help='distance metric',
    type=str,
    default='sym'
    )
parser.add_argument(
    '-s',
    '--simulations',
    dest='sims',
    help='number of simulations to do',
    type=int,
    default=1
    )
args = vars(parser.parse_args())
nsims = args['sims']
dist = args['dist']
### END ARGPARSE BIT

# Essential files
base_file = args['base']
helper = os.environ['DARWINHELPER']
output_file = args['out']
tmpdir = os.environ['TEMPORARY_DIRECTORY']

for f in [base_file, helper]:
    if not os.path.isfile(f):
        print 'Had a problem finding this file: {0}'.format(f)
        sys.exit(0)

if os.path.isfile(output_file):
    sys.exit(0)
### END FILE CHECKS

overall = []

base = cPickle.load(file(base_file))
for i in range(nsims):
    r = SequenceCollection(records=base.get_randomised_alignments(), helper=helper, tmpdir=tmpdir)
    r.put_trees_parallel(program='phyml',model='JC69',datatype='nt',ncat=1,tmpdir=tmpdir)
    r.put_partitions(metrics=[dist],linkages=['ward'],nclasses=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
    r.put_clusters()
    r.put_cluster_trees_parallel(program='phyml',model='JC69',datatype='nt',ncat=1,tmpdir=tmpdir)
    results = r.get_clusters()
    scores = []
    for j in range(1,21):
        try:
            scores.append( results[ (dist, 'ward', j)].score )
        except KeyError:
            continue
    overall.append(scores)
    
writer=open(output_file,'w')
for l in overall:
	for score in l:
		writer.write('{0}\t'.format(score))
	writer.write('\n')

writer.flush()
writer.close()
