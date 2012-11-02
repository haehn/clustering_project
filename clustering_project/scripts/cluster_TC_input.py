#!/usr/bin/env python

import argparse
import glob
import os
import re
import sys

from sequence_collection import SequenceCollection
from sequence_record import TCSeqRec
from Tree import tree


# some definitions

def fpath(s):
    """
    Helper function used when passing filepath arguments with argparse module.
    Trims all '/' characters from the end of the path string.
    """

    while s.endswith('/'):
        s = s[:-1]
    return s


sort_key = lambda item: tuple((int(num) if num else alpha) for (num,
                              alpha) in re.findall(r'(\d+)|(\D+)',
                              item))

# Parse command-line arguments

desc = \
    '''
Script to run tree-clustering analysis on TreeCollection data -

Requires:
Directory to work on, containing:
  -  one distance-variance matrix and one genome map per gene,
     in the format (gene_name).dv and (gene_name).gm
  -  one or more guide tree files, in newick format
     (the script loads all *.nwk files as guide trees, and uses them all)
  -  one set of taxon labels, Labels.txt
'''
name_of_this_program = sys.argv[0]
if '/' in name_of_this_program:
    i = name_of_this_program.rindex('/')
    name_of_this_program = name_of_this_program[i + 1:]

home_dir = os.environ['HOME']

parser = argparse.ArgumentParser(prog=name_of_this_program,
                                 description=desc)

parser.add_argument('-d', '--directory', help='Working directory',
                    type=fpath, default=home_dir)

args = vars(parser.parse_args())
working_dir = args['directory']

if not os.path.isdir(working_dir):
    print 'Can\'t find directory: {0}'.format(working_dir)
    sys.exit()

print working_dir

# Do work...

dv_files = glob.glob('{0}/*.dv'.format(working_dir))
gm_files = glob.glob('{0}/*.gm'.format(working_dir))
labels_file = glob.glob('{0}/Labels.txt'.format(working_dir))
tree_files = glob.glob('{0}/*.nwk'.format(working_dir))

for (i, fileset) in enumerate((dv_files, gm_files, labels_file,
                              tree_files)):
    if len(fileset) == 0:
        if i == 0:
            print 'Can\'t find *.dv files in {0}'.format(working_dir)
        elif i == 1:
            print 'Can\'t find *.gm files in {0}'.format(working_dir)
        elif i == 2:
            print 'Can\'t find Labels.txt in {0}'.format(working_dir)
        elif i == 3:
            print 'Can\'t find *.nwk files in {0}'.format(working_dir)
        sys.exit()

dv_files.sort(key=sort_key)
gm_files.sort(key=sort_key)
tree_files.sort(key=sort_key)
