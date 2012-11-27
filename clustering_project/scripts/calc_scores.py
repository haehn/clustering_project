#!/usr/bin/env python

from sequence_collection import SequenceCollection
from partition import Partition
import os
import sys
import argparse


def fpath(s):
    """
    Helper function used when passing filepath arguments with argparse module.
    Trims all '/' characters from the end of the path string.
    """

    while s.endswith('/'):
        s = s[:-1]
    return s


desc = 'Calculates likelihood of clusters of alignments'
parser = argparse.ArgumentParser(prog=sys.argv[0], description=desc,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-d', '--indir', help='Input directory', type=fpath)
parser.add_argument('-f', '--format', help='Alignment file format',
                    type=str, default='phylip')
parser.add_argument('-s', '--sequencedatatype',
                    help="'dna' or 'protein'", type=str, default='dna')
parser.add_argument('-n', '--nclasses', help='Number of classes',
                    type=int, default=4)
parser.add_argument('-o', '--outfile', help='Where to write results',
                    type=fpath, default='scores.txt')

args = vars(parser.parse_args())
indir = args['indir']
seqdir = '{0}/dna_alignments'.format(indir)
format = args['format']
outf = '/'.join((indir,args['outfile']))
nclasses = args['nclasses']
datatype = args['sequencedatatype']
try:
    TMPDIR = os.environ['TEMPORARY_DIRECTORY']
except:
    TMPDIR = '/tmp'
HELPER = os.environ['DARWINHELPER']
GTP_PATH = os.environ['GTP_PATH']
dists = ['rf', 'euc', 'geo']
methods = [
    'single',
    'complete',
    'average',
    'ward',
    'kmedoids',
    'spectral',
    'MDS',
    ]
calc_varinf = False

if os.path.isfile(outf):
    sys.exit(1)

sc = SequenceCollection(
    seqdir,
    tmpdir=TMPDIR,
    gtp_path=GTP_PATH,
    helper=HELPER,
    file_format=format,
    datatype=datatype,
    parallel_load=False,
    get_distances=False,
    )

sc.put_trees(program='bionj', model='GTR', tmpdir=TMPDIR, ncat=4,
             datatype='nt')

# sc.put_distance_matrices(['rf'])
# print sc.get_distance_matrices()['rf']
# sys.exit()

sc.put_partitions('rf', methods, nclasses, recalculate=True)
sc.put_partitions('euc', methods, nclasses, recalculate=True)
sc.put_partitions('geo', methods, nclasses, recalculate=True)

if os.path.isfile('{0}/treedistances.txt'.format(indir)):
    calc_varinf = True
    with open('{0}/treedistances.txt'.format(indir)) as truthf:
        true_clustering = \
            tuple(eval(truthf.readline().strip().split('\t')[1]))
    sc.clusters_to_partitions[('true', '')] = true_clustering
    sc.partitions[true_clustering] = Partition(true_clustering)

sc.concatenate_records()

sc.put_cluster_trees(program='bionj', datatype='nt', model='GTR',
                     tmpdir=TMPDIR, ncat=4)

d = sc.clusters_to_partitions

print 'writing to', outf
with open(outf, 'w') as writer:
    for k in sorted(d):
        partition_object = sc.partitions[d[k]]
        score = partition_object.score
        if calc_varinf:
            varinf = partition_object.variation_of_information(d[k],
                    true_clustering)
        else: varinf = ''

        head = ' '.join(k[:2])
        writer.write(','.join([head, str(score), str(varinf)])+'\n')
