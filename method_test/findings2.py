#!/usr/bin/env python

import glob
import numpy as np
from pylab import *
import argparse
import re
import cPickle
import sys
from collections import defaultdict
from clustering import Clustering

def fpath(s):
    """
    Helper function used when passing filepath arguments with argparse module.
        Trims all '/' characters from the end of the path string.
    """
    while s.endswith('/'):
        s = s[:-1]
    return s

getname = lambda x: x[x.rindex('/')+1:x.rindex('.')]
sort_key = lambda item: tuple((int(num) if num else alpha) for (num,alpha) in re.findall(r'(\d+)|(\D+)', item))

parser = argparse.ArgumentParser(prog='findings.py')
parser.add_argument('-d', '--directory', help='directory to search for results', type=fpath, default='.')
args = vars(parser.parse_args())
indir = args['directory']

simdirs = sorted(glob.glob('{0}/sim*'.format(indir)), key=sort_key)
phymlresults = defaultdict(list)
bionjresults = defaultdict(list)

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
 ('sym', 'ward', 4)]

for simdir in simdirs:
    print 'Loading pickles from {0}'.format(simdir)
    phyml_sc = cPickle.load(open('{0}/phyml_sc_result.pickle'.format(simdir)))
    bionj_sc = cPickle.load(open('{0}/bionj_sc_result.pickle'.format(simdir)))
    truepartition = phyml_sc.clustering.partitions['true']
    assert truepartition == bionj_sc.clustering.partitions['true']
    for k in keys:
        mlpartition = phyml_sc.clustering.partitions[k]
        njpartition = bionj_sc.clustering.partitions[k]
        mlvarinf = Clustering().variation_of_information(truepartition, mlpartition)
        njvarinf = Clustering().variation_of_information(truepartition, njpartition)
        mlscore = phyml_sc.get_clusters()[k].score
        njscore = bionj_sc.get_clusters()[k].score
        phymlresults[k].append( (mlscore, mlvarinf) )
        bionjresults[k].append( (njscore, njvarinf) )

with open('{0}/phymlresults.txt'.format(indir),'w') as file:
    for k in keys:
        mean_score = np.mean([x for x,y in phymlresults[k]])
        mean_score_stderr = np.std([x for x,y in phymlresults[k]],ddof=1)/np.sqrt(len(bionjresults[k]))
        mean_varinf = np.mean([y for x,y in phymlresults[k]])
        mean_varinf_stderr = np.std([y for x,y in phymlresults[k]],ddof=1)/np.sqrt(len(bionjresults[k]))
        file.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(' '.join(k[:2]), mean_score, mean_score_stderr, mean_varinf, mean_varinf_stderr))

with open('{0}/bionjresults.txt'.format(indir),'w') as file:
    for k in keys:
        mean_score = np.mean([x for x,y in bionjresults[k]])
        mean_score_stderr = np.std([x for x,y in bionjresults[k]],ddof=1)/np.sqrt(len(bionjresults[k]))
        mean_varinf = np.mean([y for x,y in bionjresults[k]])
        mean_varinf_stderr = np.std([y for x,y in bionjresults[k]],ddof=1)/np.sqrt(len(bionjresults[k]))
        file.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(' '.join(k[:2]), mean_score, mean_score_stderr, mean_varinf, mean_varinf_stderr))
