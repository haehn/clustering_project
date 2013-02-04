#!/usr/bin/env python
print 'going'
################################################################################
# methodTest.py
# -------------

# Script to run clustering procedure on input data using all
# implemented clustering methods and all implemented distance metrics
# to decide which one is best (SEE ALSO calc_scores.py)

# command-line args:
#    -i input directory
#    -o output scores file
#    -f file format (['fasta', 'phylip'])
#    -t tree-building method to use (['phyml', 'bionj', 'treecollection'])
#    -d datatype (['protein', 'dna'])
#    -n nclasses
#    -l level
#    -s simulation rep
#    -g number of guidetrees to use (treecollection)
################################################################################

from sequence_collection import SequenceCollection
from partition import Partition
import sys
import os
import argparse


def topdir(s):
    if '/' in s:
        return s[s.rindex('/') + 1:]
    return s

prog = topdir(sys.argv[0])
desc = \
    '''\
Script to run clustering procedure on input data using all 
implemented clustering methods and all implemented distance metrics
to decide which one is best'''

tchoices = ['phyml', 'bionj', 'treecollection']
dchoices = ['protein', 'dna']
fchoices = ['fasta', 'phylip']

parser = argparse.ArgumentParser(prog=prog, description=desc)
parser.add_argument('-i', dest='indir', type=str)
parser.add_argument('-o', dest='outf', type=str)
parser.add_argument('-f', dest='format', choices=fchoices, type=str)
parser.add_argument('-t', dest='treeprog', choices=tchoices, type=str)
parser.add_argument('-d', dest='datatype', choices=dchoices, type=str)
parser.add_argument('-n', dest='nclasses', type=int, nargs='*',
                    default=4)
parser.add_argument('-l', dest='level', type=int, required=True)
parser.add_argument('-s', dest='sim', type=int, required=True)
parser.add_argument('-g', dest='guide', type=int, default=10)
args = parser.parse_args()

indir = args.indir.rstrip('/')
outf = args.outf.rstrip('/')
treeprog = args.treeprog
datatype = args.datatype
file_format = args.format
nclasses = args.nclasses
level = args.level
rep = args.sim
nguide = args.guide

# print indir, treeprog, datatype, file_format, nclasses
# sys.exit()

if not 'GTP_PATH' in os.environ:
    gtp_path = '/'.join([
        '',
        'nfs',
        'research2',
        'goldman',
        'kevin',
        'git',
        'kevin',
        'clustering_project',
        'class_files',
        ])
else:
    gtp_path = os.environ['GTP_PATH']

if not 'DARWINHELPER' in os.environ:
    helper = '/'.join([
        '',
        'nfs',
        'research2',
        'goldman',
        'kevin',
        'git',
        'kevin',
        'clustering_project',
        'class_files',
        'DV_wrapper.drw',
        ])
else:
    helper = os.environ['DARWINHELPER']

try:
    TMPDIR = os.environ['TEMPORARY_DIRECTORY']
except:
    TMPDIR = '/tmp'

### MAIN
if treeprog == 'treecollection':
    get_distances = True
else:
    get_distances = False
sc = SequenceCollection(indir, file_format=file_format,
                        gtp_path=gtp_path, datatype=datatype,
                        get_distances=get_distances, tmpdir=TMPDIR, helper=helper)

sc.put_trees(program=treeprog)
sc.put_partitions(['geo', 'euc', 'rf'], [
    'average',
    'complete',
    'kmedoids',
    'MDS',
    'single',
    'spectral00',
    'spectral01',
    'spectral10',
    'spectral11',
    'ward',
    ], nclasses, recalculate=True)

if os.path.isfile('{0}/../treedistances.txt'.format(indir)):
    calc_varinf = True
    with open('{0}/../treedistances.txt'.format(indir)) as truthf:
        true_clustering = \
            tuple(eval(truthf.readline().strip().split('\t')[1]))
    sc.clusters_to_partitions[('true', 'na')] = true_clustering
    sc.partitions[true_clustering] = Partition(true_clustering)

sc.concatenate_records()
sc.put_cluster_trees(program=treeprog, tmpdir=TMPDIR, max_guide_trees=nguide)

d = sc.clusters_to_partitions
true_score = sc.partitions[sc.clusters_to_partitions[('true','na')]].score
print 'writing to', outf
with open(outf, 'w') as writer:
    for k in sorted(d):
        if k[0] == 'true': continue
        partition_object = sc.partitions[d[k]]
        score = partition_object.score
        if calc_varinf:
            varinf = partition_object.variation_of_information(d[k],
                    true_clustering)
        else:
            varinf = ''
        
        head = ','.join(k[:2])
        writer.write(','.join([head, str(score), str(round(score-true_score, 10)), str(round(varinf,10)), str(level), str(rep)]) + '\n')
