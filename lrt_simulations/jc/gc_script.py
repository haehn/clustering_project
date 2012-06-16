#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import cPickle
import os
import sys

# Command-line arguments handled by argparse module
parser = argparse.ArgumentParser(prog='gc_script.py',
                description='Run Goldman-Cox procedure on SequenceCollection'
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
    help='output pickle file',
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
    '-m',
    '--merged',
    dest='merged',
    help='Set this flag to use merged clusters',
    action='store_true',
    default=False
    )
parser.add_argument(
    '-s',
    '--simulations',
    dest='sims',
    help='number of simulations to do',
    type=int,
    default=1
    )
parser.add_argument(
    '-c',
    '--classes',
    dest='classes',
    help='number of classes to find',
    type=int,
    default=1
    )

args = vars(parser.parse_args())
nsims = args['sims']
nclasses = args['classes']
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
### END FILE CHECKS

# main
results = []
base = cPickle.load(file(base_file))
for i in range(nsims):
    simulated = base.simulate_from_result((dist, 'ward', nclasses), helper=helper, tmpdir=tmpdir)
    simulated.put_trees_parallel(program='phyml', model='JC69', ncat=1, datatype='nt', tmpdir=tmpdir)
    simulated.put_partitions(dist,'ward',[nclasses,nclasses+1])
    simulated.put_clusters()
    simulated.put_cluster_trees_parallel(program='phyml', model='JC69', ncat=1, datatype='nt', tmpdir=tmpdir)
    r1 = simulated.get_clusters()[(dist,'ward',nclasses)]
    r2 = simulated.get_clusters()[(dist,'ward',nclasses+1)]
    diff = r2.score-r1.score
    results.append(diff)
writer = open(output_file,'w')
for result in results:
    writer.write('{0}\n'.format(result))
writer.flush()
writer.close()
print results
