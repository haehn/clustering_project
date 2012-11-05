#!/usr/bin/env python

import argparse
import glob
import os
import re
import sys

from sequence_collection import SequenceCollection
from sequence_record import TCSeqRec
from tree import Tree


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

get_name = lambda filename: filename[filename.rindex('/')
    + 1:filename.rindex('.')]


def get_best_TC_tree(
    dv_file,
    gm_file,
    label_file,
    tree_files,
    name='unnamed_tree',
    ):
    """
    Given a distance-variance file, a genome-map file, a label file
    and a number of tree files
    """

    if not isinstance(tree_files, list):
        return Tree.new_treecollection_tree(dv_file, gm_file,
                label_file, tree_files, name)
    best_score = float('inf')
    best_tree = None
    starts = len(tree_files)

    for i in range(starts):
        guide = tree_files[i]
        current_tree = Tree.new_treecollection_tree(dv_file, gm_file,
                label_file, guide, name)
        if current_tree.score < best_score:
            best_score = current_tree.score
            best_tree = current_tree

    return best_tree

##############################
# Parse command-line arguments
##############################

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

print 'Working directory...', working_dir

###################
# Do file-gathering
###################

dv_files = glob.glob('{0}/*.dv'.format(working_dir))
gm_files = glob.glob('{0}/*.gm'.format(working_dir))
labels_file = glob.glob('{0}/Labels.txt'.format(working_dir))
tree_files = glob.glob('{0}/*.nwk'.format(working_dir))

# Check for presence of files
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

# Check number of dv_files == number of gm_files
if len(dv_files) != len(gm_files):
    print 'Number of distance-variance files doesn\'t match'
    print 'number of genome map files'
    sys.exit()

dv_files.sort(key=sort_key)
gm_files.sort(key=sort_key)
tree_files.sort(key=sort_key)
labels_file = labels_file[0]

#######################
# Make TCSeqRec objects
#######################

number_of_genes = len(dv_files)
records = []
for i in range(number_of_genes):
    print i
    dv = dv_files[i]
    gm = gm_files[i]
    # read file contents
    with open(dv) as dv_file_reader:
        dv_matrix = dv_file_reader.read()
    with open(labels_file) as labels_file_reader:
        labels = labels_file_reader.read()
    name = get_name(dv)
    tree = get_best_TC_tree(dv, gm, labels_file, tree_files, name)

    record = TCSeqRec()
    record.dv.append((dv_matrix, labels))
    record.tree = tree
    record.name = name
    records.append(record)





