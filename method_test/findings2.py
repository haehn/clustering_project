#!/usr/bin/env python

import argparse
import glob
from collections import defaultdict

def fpath(s):
    """
    Helper function used when passing filepath arguments with argparse module.
        Trims all '/' characters from the end of the path string.
    """
    while s.endswith('/'):
        s = s[:-1]
    return s

parser = argparse.ArgumentParser(prog='findings.py')
parser.add_argument('-d', '--directory', help='directory to search for results', type=fpath, default='.')
args = vars(parser.parse_args())
basedir = args['directory']

phymldic = defaultdict(list)
bionjdic = defaultdict(list)
phymlvidic = defaultdict(list)
bionjvidic = defaultdict(list)

simdirs = glob.glob('{0}/sim*'.format(basedir))

print 'Reading results...'
for simdir in simdirs:
    with open('{0}/phyml_results.txt'.format(simdir)) as phyml_results:
        for line in phyml_results:
            name, score, varinf = line.rstrip().split('\t')
            phymldic[name].append(score)
            phymlvidic[name].append(varinf)

    with open('{0}/bionj_results.txt'.format(simdir)) as bionj_results:
        for line in bionj_results:
            name, score, varinf = line.rstrip().split('\t')
            bionjdic[name].append(score)
            bionjvidic[name].append(varinf)
print 'Done.'

keys = ['true',
 'euc MDS',
 'euc average',
 'euc complete',
 'euc kmedoids',
 'euc single',
 'euc spectral',
 'euc ward',
 'geo MDS',
 'geo average',
 'geo complete',
 'geo kmedoids',
 'geo single',
 'geo spectral',
 'geo ward',
 'sym MDS',
 'sym average',
 'sym complete',
 'sym kmedoids',
 'sym single',
 'sym spectral',
 'sym ward'
 ]

print 'Writing summaries...'

with open('{0}/phyml_lnl_results.txt'.format(basedir),'w') as phylnlout:
    for k in keys:
        phylnlout.write('{0}\t'.format(k))
        phylnlout.write(' '.join(phymldic[k]))
        phylnlout.write('\n')

with open('{0}/phyml_vi_results.txt'.format(basedir),'w') as phyviout:
    for k in keys:
        phyviout.write('{0}\t'.format(k))
        phyviout.write(' '.join(phymlvidic[k]))
        phyviout.write('\n')

with open('{0}/bionj_lnl_results.txt'.format(basedir),'w') as njlnlout:
    for k in keys:
        njlnlout.write('{0}\t'.format(k))
        njlnlout.write(' '.join(bionjdic[k]))
        njlnlout.write('\n')

with open('{0}/bionj_vi_results.txt'.format(basedir),'w') as njviout:
    for k in keys:
        njviout.write('{0}\t'.format(k))
        njviout.write(' '.join(bionjvidic[k]))
        njviout.write('\n')   
print 'Done.'

