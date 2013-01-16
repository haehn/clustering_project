#!/usr/bin/env python

################################################################################
# Read in sequence alignments and prepared trees (optionally can delete these)
# Do clustering and concatenate records
# Dump concatenated records to disk, and pickle record
################################################################################

################################################################################
# Commandline args:
#       -i      = input directory containing alignments and trees
#       -o      = output filename
#       -d      = distance method to use
#       -c      = clustering method to use
#       -min    = minimum number of clusters
#       -max    = maximum number of clusters
#       -t      = temp directory
#       -del    = delete alignments and trees
#       -data   = datatype (protein or dna)
#       -score  = evaluate cluster trees in this script, not later in LSF
################################################################################

import argparse
import glob
import os
import re
import sys
from errors import directorycheck_and_quit, directorycheck_and_make

progname = re.compile('[A-Za-z0-9.-_]+').search(sys.argv[0]).group()
desc = '\n'.join([progname.upper(),
                 'Read in sequence alignments and prepared trees',
                 'Do clustering and concatenate records',
                 'Dump concatenated records to disk, and pickle record'
                 ])
input_help = 'Path to the input directory'
output_help = 'Name of output file'
tmpdir_help = 'Directory used for temp files'
min_help = 'Minimum number of clusters to find'
max_help = '\n'.join(['Maximum number of clusters to find. All',
                     'intermediate numbers between \'min\' and',
                     '\'max\' will be found too.'])
distance_help = 'Distance metric(s) to use'
method_help = 'Clustering method(s) to use'
delete_help = 'Deletes input files after reading (data is saved in '
datatype_help = 'Datatype = \'protein\' || \'dna\''
valid_datatypes = ['protein', 'dna']
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
parser = argparse.ArgumentParser(prog=progname, description=desc,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', dest='input', help=input_help, type=str)
parser.add_argument('-o', dest='output', help=output_help, type=str)
parser.add_argument('-t', dest='tmpdir', help=tmpdir_help, type=str)
parser.add_argument('-min', dest='min_clusters', help=min_help,
                    type=int, default=2)
parser.add_argument('-max', dest='max_clusters', help=max_help,
                    type=int, default=20)
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
parser.add_argument('-del', dest='delete', action='store_true',
                    help=delete_help)
parser.add_argument('-data', dest='datatype', choices=valid_datatypes,
                    type=str, help=datatype_help)
parser.add_argument('-score', dest='score', action='store_true')

args = vars(parser.parse_args())
input_dir = args['input'].rstrip('/')
output = args['output']
tmpdir = args['tmpdir'].rstrip('/')
min_clusters = args['min_clusters']
max_clusters = args['max_clusters']
distance = args['distance']
method = args['cluster_method']
delete = args['delete']
datatype = args['datatype']
score = args['score']

directorycheck_and_quit(input_dir)
directorycheck_and_make(tmpdir)

gtp_path = os.environ['GTP_PATH']
helper = os.environ['DARWINHELPER']

from sequence_collection import SequenceCollection

sc = SequenceCollection(
    input_dir,
    file_format='phylip',
    datatype=datatype,
    helper=helper,
    gtp_path=gtp_path,
    tmpdir=tmpdir,
    overwrite=True,
    )

sc.load_phyml_results(input_dir, program=None)
sc.quality_scores = {}
for dist in distance:
    (_, qs) = sc.autotune(dist, max_groups=max_clusters,
                          min_groups=min_clusters)
    sc.quality_scores[dist] = qs
cluster_range = range(min_clusters, max_clusters+1)
if min_clusters > 1:
    cluster_range.insert(0, 1)
sc.put_partitions(distance, method, cluster_range)
sc.concatenate_records()
if score:
       sc.put_cluster_trees(program='bionj', optimise='r', ncat=4 )

try:
    sc.gzip(output)
except:
    print 'Couldn\'t save pickle'
    raise
    sys.exit()

if delete:
    for f in glob.glob('{0}/*.phy'.format(input_dir)) \
        + glob.glob('{0}/*_phyml_*'.format(input_dir)):
        os.remove(f)
