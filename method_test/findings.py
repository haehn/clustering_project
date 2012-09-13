#!/usr/bin/python

print 'Importing modules...'
# IMPORTS
import argparse
import cPickle
import glob
import re
import sys
import os
from sequence_collection import SequenceCollection
# END IMPORTS
print '(Done).'

# DEFINITIONS
def fpath(s):
    """
    Helper function used when passing filepath arguments with argparse module.
        Trims all '/' characters from the end of the path string.
    """
    while s.endswith('/'):
        s = s[:-1]
    return s

def load_records(directory, extension='*.pickle'):
    return [cPickle.load(open(f)) for f in
            glob.iglob('{0}/{1}'.format(directory, extension))]

def make_target_dict(target_list):
    return dict( (rec.name, rec) for rec in target_list)

def rebuild_partitions(query_list, target_dict, metric, method):
    """
    The query list is the individual alignments record list,
    The target list is the clustered alignments record list
    """

    partition_list = [0] * len(query_list)
    partition_list_check_copy = [0] * len(query_list)
    header = query_list[0].headers[0]
    header_check_copy = query_list[0].headers[-1]
    key = metric[:3] + method[:3] + '4_'

    for i in range(1,5):
        for j, rec in enumerate(query_list):
            query_seq = rec.mapping[header]
            query_seq_check_copy = rec.mapping[header_check_copy]
            if query_seq in target_dict[key+str(i)].mapping[header]:
                partition_list[j] = i
            if query_seq_check_copy in target_dict[key+str(i)].mapping[header_check_copy]:
                partition_list_check_copy[j] = i

    try:
        assert partition_list == partition_list_check_copy
    except AssertionError:
        print 'Can\'t get partitions to match'
        return -1
    return partition_list

sort_key = lambda item: tuple((int(num) if num else alpha) for (num,alpha) in re.findall(r'(\d+)|(\D+)', item.name))
# END DEFINITIONS

# COMMAND-LINE PARSING
parser = argparse.ArgumentParser(prog='findings.py')
parser.add_argument('-d', '--directory', help='directory to search for results', type=fpath, default='.')
parser.add_argument('-r', '--record_dir', help='directory in which the initial tree inference records are pickled', type=fpath, default='dna_alignments')
parser.add_argument('-p', '--program', help='Which program was used to do tree inference (phyml or bionj)', type=str, default='phyml')
parser.add_argument('-gtp', help='path to gtp.jar', type=fpath, default='homes/kgori/research/')
args = vars(parser.parse_args())
simdir = args['directory']
record_dir = args['record_dir']
program = args['program']

try: index=os.environ['LSB_JOBINDEX']
except KeyError: index = None

if index: simdir += index

try: tmpdir = os.environ['TEMPORARY_DIRECTORY']
except KeyError: tmpdir = '/tmp'
# END COMMAND-LINE PARSING

# MAIN

print 'Working directory = ',simdir
print 'Records located in',record_dir
print 'Program = ', program
print 'Temp directory =', tmpdir
print 'Reading records...'

# some initialisation
if program == 'bionj':
    records = sorted(load_records('/'.join((simdir, record_dir)), '*.nj.pickle'), key=sort_key)
    cluster_records = sorted(load_records('/'.join((simdir, 'bionj_clustering'))), key=sort_key)
else:
    records = sorted(load_records('/'.join((simdir, record_dir)), '*.ml.pickle'), key=sort_key)
    cluster_records = sorted(load_records('/'.join((simdir, 'phyml_clustering'))), key=sort_key)

cluster_dic = make_target_dict(cluster_records)
print '(Done).'

#rebuild sequenceCollection object

sc = SequenceCollection()
sc.records = records

print 'Generating distance matrices...'
for metric in ['euc', 'sym', 'geo']:
    sc.put_distance_matrices(metric, gtp_path = '/homes/kgori/research/clustering_project/class_files', tmpdir=tmpdir)
    for method in ['single', 'complete', 'ward', 'average', 'spectral', 'MDS', 'kmedoids']:
        partition = rebuild_partitions(records, cluster_dic, metric=metric, method=method)
        sc.clustering.partitions[(metric, method, 4)] = partition

sc.clustering.partitions['true'] = [1]*15 + [2]*15 + [3]*15 + [4]*15

sc.put_clusters()
print '(Done).'

for rec in sc.get_cluster_records():
    try:
        rec.tree = cluster_dic[rec.name].tree
    except KeyError:
        new_key = 'true_cluster{0}'.format(rec.name[-1])
        rec.tree = cluster_dic[new_key].tree

results_dic = sc.get_clusters()

keys = ['true',
        ('euc', 'MDS', 4),
        ('euc', 'average', 4),
        ('euc', 'complete', 4),
        ('euc', 'kmedoids', 4),
        ('euc', 'single', 4),
        ('euc', 'spectral', 4),
        ('euc', 'ward', 4),
        ('geo', 'MDS', 4),
        ('geo', 'average', 4),
        ('geo', 'complete', 4),
        ('geo', 'kmedoids', 4),
        ('geo', 'single', 4),
        ('geo', 'spectral', 4),
        ('geo', 'ward', 4),
        ('sym', 'MDS', 4),
        ('sym', 'average', 4),
        ('sym', 'complete', 4),
        ('sym', 'kmedoids', 4),
        ('sym', 'single', 4),
        ('sym', 'spectral', 4),
        ('sym', 'ward', 4)
        ]

print 'Writing results to {0}/{1}_results.txt'.format(simdir, program)
with open('{0}/{1}_results.txt'.format(simdir, program), 'w') as output_file:
    for k in keys:
        result = results_dic[k]
        result.update_score()
        score = result.score
        varinf = sc.clustering.variation_of_information(sc.clustering.partitions[k], sc.clustering.partitions['true'])
        if type(k) is tuple:
            output_file.write('{0}\t{1}\t{2}\n'.format(' '.join(k[:2]), score, varinf))
        else:
            output_file.write('{0}\t{1}\t{2}\n'.format(k, score, varinf))
print '(All done).'