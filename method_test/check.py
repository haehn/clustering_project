#!/usr/bin/env python

import glob
import argparse
import os


def fpath(s):
    """
    Helper function used when passing filepath arguments with argparse module.
        Trims all '/' characters from the end of the path string.
    """

    while s.endswith('/'):
        s = s[:-1]
    return s

parser = argparse.ArgumentParser(prog='check.py')
parser.add_argument('-d', '--base-directory', type=fpath)
args = vars(parser.parse_args())
basedir = args['base_directory']

simdirs = set(glob.glob('{0}/sim*'.format(basedir)))
expected_dirs = set(('{0}/sim{1}'.format(basedir, x) for x in range(1,101)))
diff = expected_dirs-simdirs
if diff: print diff

exp_results = set(('{0}{1}.txt'.format(a,b) for a in ['euc','geo','sym'] for b in ['single','complete','average','ward','MDS','kmedoids','spectral']))
for directory in simdirs:
	results = glob.glob('{0}/euc*'.format(directory)) \
        + glob.glob('{0}/geo*'.format(directory)) \
        + glob.glob('{0}/sym*'.format(directory))
	results = set([x[x.rindex('/')+1:] for x in results])
	diff2 = exp_results - results
	if diff2:
		print directory, diff2

