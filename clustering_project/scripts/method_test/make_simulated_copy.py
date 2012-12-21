#!/usr/bin/env python

################################################################################
# Read in a SequenceCollection from disk and dump a simulated copy
################################################################################

################################################################################
# Commandline args:
#       -i      = input filename
#       -o      = output filename
#       -d      = distance
#       -c      = clustering method
#       -m      = tree-building method (phyml || bionj || bionj+)
#       -ind    = simulation number (index)
#       -min    = minimum number of clusters
#       -max    = maximum number of clusters
#       -t      = temp dir
################################################################################

import argparse
import os 
import re
import sys
from errors import filecheck_and_quit, directorycheck_and_make

progname = re.compile('[A-Za-z0-9.-_]+').search(sys.argv[0]).group()
desc = \
    'Read in a SequenceCollection from disk and dump a simulated copy'
input_help = 'Filepath+name of gzipped SequenceCollection object'
output_help = 'Output path'
min_help = 'Minimum number of clusters to find'
max_help = '\n'.join(['Maximum number of clusters to find. All',
                     'intermediate numbers between \'min\' and',
                     '\'max\' will be found too.'])
distance_help = 'Distance metric(s) to use'
method_help = 'Clustering method(s) to use'

valid_dists = ['geo', 'euc', 'rf']
valid_methods = [
    'spectral',
    'MDS',
    'ward',
    'kmedoids',
    'single',
    'complete',
    'average',
    ]
valid_phylo = ['phyml', 'bionj', 'bionj+']

parser = argparse.ArgumentParser(prog=progname, description=desc)
parser.add_argument('-i', dest='input_file', help=input_help, type=str)
parser.add_argument('-o', dest='output_dir', help=output_help, type=str)
parser.add_argument('-t', dest='tmpdir', type=str, default=None)
parser.add_argument('-min', dest='min_clusters', help=min_help,
                    type=int, default=2)
parser.add_argument('-max', dest='max_clusters', help=max_help,
                    type=int, default=20)
parser.add_argument('-ind', dest='ind', type=int, default=1)
parser.add_argument(
    '-m',
    dest='tree_method',
    help=method_help,
    type=str,
    default='bionj+',
    choices=valid_phylo,
    )
parser.add_argument(
    '-d',
    dest='distance',
    help=distance_help,
    type=str,
    choices=valid_dists,
    nargs='+',
    default='geo',
    )
parser.add_argument(
    '-c',
    dest='cluster_method',
    help=method_help,
    type=str,
    choices=valid_methods,
    nargs='+',
    default='spectral',
    )

args         = parser.parse_args()
input_file   = args.input_file
output_dir   = args.output_dir
tmpdir       = args.tmpdir
min_clusters = args.min_clusters
max_clusters = args.max_clusters
distance     = args.distance
method       = args.cluster_method
ind          = args.ind
tree_method  = args.tree_method

filecheck_and_quit(input_file)
directorycheck_and_make(output_dir)
directorycheck_and_make(tmpdir)

################################################################################
# Main 
################################################################################

from sequence_collection import SequenceCollection

sc = SequenceCollection.gunzip(input_file)
for c in range(min_clusters, max_clusters):
    try:
        assert (distance, method, c) in sc.clusters_to_partitions
    except AssertionError:
        print c
        sys.exit()

for c in range(min_clusters, max_clusters):
    output_subdir = '{0}/s{1}_cl{2}'.format(output_dir, ind, c)
    directorycheck_and_make(output_subdir, verbose=False)
    print 'Simulating {0}-part clustering from {1}'.format(c,
            input_file)
    sc.simulate_from_result((distance, method, c),
                            output_dir=output_subdir, name='sim',
                            tmpdir=tmpdir)
    if os.environ['SCRIPTS']:
        script_dir = os.environ['SCRIPTS']
        os.system('{0}/lsf_phyml.py -i {1} -d {2} -m {3}'.format(script_dir,
                  output_subdir, sc.datatype, tree_method))

################################################################################
# NEXT STEPS: 
    # $SCRIPTS/read_alignments_and_cluster.py -i output_subdir -score [... -del]
    # $SCRIPTS/extract_scores.py [...]
################################################################################