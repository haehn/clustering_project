#!/usr/bin/env python

import glob
import numpy as np
from pylab import *
import argparse
import re
import cPickle
import sys
import os

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
simdir = args['directory']

try: index=os.environ['LSB_JOBINDEX']
except: index = None

if index: simdir += index

print 'Loading pickles from {0}'.format(simdir)
phyml_sc = cPickle.load(open('{0}/phyml_sc.pickle'.format(simdir)))
bionj_sc = cPickle.load(open('{0}/bionj_sc.pickle'.format(simdir)))
print 'Making dictionaries...'
bionjdic = dict([(getname(x),cPickle.load(open(x))) for x in glob.glob('{0}/bionj_clustering/*.pickle'.format(simdir))])
phymldic = dict([(getname(x),cPickle.load(open(x))) for x in glob.glob('{0}/phyml_clustering/*.pickle'.format(simdir))])
print '  filling in results for bionj'
for njrec in bionj_sc.get_cluster_records():
    njrec.tree = bionjdic[njrec.name].tree
print '  filling in results for phyml'
for mlrec in phyml_sc.get_cluster_records():
    mlrec.tree = phymldic[mlrec.name].tree
njresults = bionj_sc.get_clusters()
mlresults = phyml_sc.get_clusters()
print 'Updating scores'
for njk in njresults:
    njresults[njk].update_score()
for mlk in mlresults:
    mlresults[mlk].update_score()
print 'Writing results'
cPickle.dump(phyml_sc, open('{0}/phyml_sc_result.pickle'.format(simdir), 'w'))
cPickle.dump(bionj_sc, open('{0}/bionj_sc_result.pickle'.format(simdir), 'w'))

print '(done).'
