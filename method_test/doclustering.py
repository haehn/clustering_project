#!/usr/bin/env python

import argparse
from sequence_collection import SequenceCollection
import cPickle
import os
import glob
import re
import sys
import shutil

"""
Usage: bsub -o /dev/null -e error.%J.%I -J "jobname[1-x]" bash <path>/tempdir_wrapper.sh python doclustering.py -d <path>/nni/level1/sim_

Where the path endswith an underscore above, JobArray will fill in the number of the sim from 1 to x

"""

def fpath(s):
    """
    Helper function used when passing filepath arguments with argparse module.
        Trims all '/' characters from the end of the path string.
    """
    while s.endswith('/'):
        s = s[:-1]
    return s

sort_key = lambda item: tuple((int(num) if num else alpha) for (num,alpha) in re.findall(r'(\d+)|(\D+)', item))

parser = argparse.ArgumentParser(prog='doclustering.py')
parser.add_argument('-d', '--directory', help='input directory', type=fpath, default='.')
args = vars(parser.parse_args())
index = os.environ['LSB_JOBINDEX']
indir = os.path.abspath(args['directory'])
if index and index != '0':
    indir += index
print 'Working on {0}'.format(indir)
dnadir = '{0}/dna_alignments'.format(indir)
try:
    tmpdir = os.environ['TEMPORARY_DIRECTORY']
except:
    tmpdir = '/tmp'

if os.path.isdir('{0}/bionj_clustering'.format(indir)):
    shutil.rmtree('{0}/bionj_clustering'.format(indir))
if os.path.isdir('{0}/phyml_clustering'.format(indir)):
    shutil.rmtree('{0}/phyml_clustering'.format(indir))

phymlpickles = [x for x in glob.glob('{0}/*.ml.pickle'.format(dnadir)) if not '.nj.' in x]
bionjpickles = [x for x in glob.glob('{0}/*.nj.pickle'.format(dnadir)) if '.nj.' in x]
phymlrecords = sorted([cPickle.load(file(x)) for x in phymlpickles], key=lambda x:sort_key(x.name))
bionjrecords = sorted([cPickle.load(file(x)) for x in bionjpickles], key=lambda x:sort_key(x.name))

true = eval(open('{0}/treedistances.txt'.format(indir)).read().split('\n')[0].split('\t')[1])

for rec in phymlrecords:
    rec.datatype='dna'
for rec in bionjrecords:
    rec.datatype='dna'

try: 
    assert len(phymlrecords) == len(bionjrecords) == 60
except: 
    print 'Missing records in {0}'.format(indir)
    sys.exit(1)

phyml_sc = SequenceCollection(records=phymlrecords, datatype='dna', helper=os.environ['DARWINHELPER'],tmpdir=tmpdir, get_distances=False)
bionj_sc = SequenceCollection(records=bionjrecords, datatype='dna', helper=os.environ['DARWINHELPER'],tmpdir=tmpdir, get_distances=False)

phyml_sc.put_partitions(['geo','euc','sym'],['single','complete','average','ward','kmedoids','MDS','spectral'], 4, gtp_path='/net/isilon7/nobackup/research/goldman/kevin/clustering_project/class_files', tmpdir=tmpdir)
bionj_sc.put_partitions(['geo','euc','sym'],['single','complete','average','ward','kmedoids','MDS','spectral'], 4, gtp_path='/net/isilon7/nobackup/research/goldman/kevin/clustering_project/class_files', tmpdir=tmpdir)
phyml_sc.clustering.partitions['true']=true
bionj_sc.clustering.partitions['true']=true

phyml_sc.put_clusters()
bionj_sc.put_clusters()

if not os.path.isdir('{0}/phyml_clustering'.format(indir)):
    os.mkdir('{0}/phyml_clustering'.format(indir))
if not os.path.isdir('{0}/bionj_clustering'.format(indir)):
    os.mkdir('{0}/bionj_clustering'.format(indir))

phyml_result = phyml_sc.get_clusters()['true']
bionj_result = bionj_sc.get_clusters()['true']

for i, rec in enumerate(phyml_result.concats, start=1):
    rec.name = 'true_cluster{0}'.format(i)
    rec.datatype = 'dna'

for i, rec in enumerate(bionj_result.concats, start=1):
    rec.name = 'true_cluster{0}'.format(i)
    rec.datatype = 'dna'

for rec in phyml_sc.get_cluster_records():
    rec.datatype='dna'
    rec.write_phylip('{0}/phyml_clustering/{1}.phy'.format(indir,rec.name))
for rec in bionj_sc.get_cluster_records():
    rec.datatype='dna'
    rec.write_phylip('{0}/bionj_clustering/{1}.phy'.format(indir,rec.name))

cPickle.dump(phyml_sc, open('{0}/phyml_sc.pickle'.format(indir),'w')) # Let's not pickle the whole SequenceCollection, it's too big
cPickle.dump(bionj_sc, open('{0}/bionj_sc.pickle'.format(indir),'w'))

