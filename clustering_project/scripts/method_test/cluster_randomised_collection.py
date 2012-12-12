#!/usr/bin/env python

################################################################################
# Reads 
################################################################################

################################################################################
# Commandline args:
#      -i = input file: a python pickle serialization of
#                       a SequenceCollection object
#                       containing a set of sequences
#
#      -o = output directory: the randomised copies of the sequences
#                       will be dumped in here, in phylip format.
#                       Names are hashed, so a translation dictionary
#                       is pickled here, along with the randomised
#                       copy of the original SequenceCollection object
#
# Outputs:
#      <files>.phy
#      hash_translation.pkl
#      scrand.pkl
################################################################################

import cPickle
import tarfile
from errors import filecheck_and_quit, directorycheck_and_raise, \
    directorycheck_and_make_recursive
import sys
import argparse
import re
import os

progname      = re.compile('[A-Za-z0-9.-_]+').search(sys.argv[0]).group()

desc          = \
    '\n'.join(['Load the pickle containing the randomised collection,',
              'load in the precomputed phyml results,',
              'cluster the trees and calculate the scores'])
input_help    = 'Path to input directory'
tmpdir_help   = 'Path to use for tmpfiles. Will be created recursively.'
min_help      = 'Minimum number of clusters to find'
max_help      = '\n'.join(['Maximum number of clusters to find. All',
                     'intermediate numbers between \'min\' and',
                     '\'max\' will be found too.'])
distance_help = ''
method_help   = ''
valid_dists   = ['geo', 'euc', 'rf']
valid_methods = ['spectral', 'MDS', 'ward', 'kmedoids', 'single', 'complete', 
                'average']

parser = argparse.ArgumentParser(prog=progname, description=desc,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-i', '--input', help=input_help, type=str)
# parser.add_argument('-p', '--pickle', help=pickle_help, type=str)
parser.add_argument('-t', '--tmpdir', help=tmpdir_help, type=str)
parser.add_argument('-min', '--min_clusters', help=min_help, type=int,
                    default=2)
parser.add_argument('-max', '--max_clusters', help=max_help, type=int,
                    default=20)
parser.add_argument('-d', '--distance', help=distance_help, type=str,
                    choices=valid_dists, default='geo')
parser.add_argument('-c', '--cluster-method', help=method_help, type=str,
                    choices=valid_methods, default='spectral')

args         = vars(parser.parse_args())
input_dir    = args['input'].rstrip('/')
tmpdir       = args['tmpdir'].rstrip('/')
min_clusters = args['min_clusters']
max_clusters = args['max_clusters']
distance     = args['distance']
method       = args['cluster_method']
pickle = '{0}/scrand.pkl'.format(input_dir)

directorycheck_and_raise(input_dir)
directorycheck_and_make_recursive(tmpdir)
filecheck_and_quit(pickle)

sc = cPickle.load(open(pickle))
sc.tmpdir = tmpdir
print 'Loading phyml results...'
sc.load_phyml_results(input_dir, use_hashname=True)
print 'Autotuning...'
sc.autotune(distance, max_groups=max_clusters)
print 'Clustering...'
sc.put_partitions(distance, method, range(min_clusters, max_clusters))
sc.concatenate_records()
sc.put_cluster_trees(program='bionj')
scores = sorted(sc.get_scores(), key=lambda x:  x[0])
print 'Scores:'
for score in scores: print score