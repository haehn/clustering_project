#!/usr/bin/env python

import argparse
import glob
import os
import re
import sys
import cPickle

from sequence_collection import SequenceCollection
from sequence_record import TCSeqRec
from tree import Tree
from partition import Partition


# some definitions

def fpath(s):
    """
    Helper function used when passing filepath arguments with argparse module.
    Trims all '/' characters from the end of the path string.
    """

    return s.rstrip('/')


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


def multiwordReplace(text, wordDic):
    """
    take a text and replace words that match a key in a dictionary with
    the associated value, return the changed text
    """

    rc = re.compile('|'.join(map(re.escape, wordDic)))

    def translate(match):
        return wordDic[match.group(0)]

    return rc.sub(translate, text)


names = dict(
    HUMAN='Human',
    MOUSE='Mouse',
    CANFA='Dog',
    RATNO='Rat',
    BOVIN='Cow',
    PANTR='Chimp',
    MONDO='Opossum',
    MACMU='RhesusMonkey',
    LOXAF='Elephant',
    RABIT='Rabbit',
    ECHTE='LesserHedgehogTenrec',
    OTOGA='NorthernGreaterGalago',
    TUPGB='NorthernTreeShrew',
    MYOLU='LittleBrownBat',
    CAVPO='GuineaPig',
    MICMU='GreyMouseLemur',
    PONAB='Orangutan',
    HORSE='Horse',
    TURTR='BottlenoseDolphin',
    PTEVA='LargeFlyingFox',
    PROCA='RockHyrax',
    DIPOR='KangarooRat',
    TARSY='PhilippineTarsier',
    GORGO='Gorilla',
    CALJA='CommonMarmoset',
    AILME='Panda',
    SARHA='TasmanianDevil',
    NOMLE='NorthernWhiteCheekedGibbon',
    )

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
parser.add_argument('-o', '--output',
                    help='Output filename (saves in -d path)',
                    type=str, default='collection')

args = vars(parser.parse_args())
working_dir = args['directory']
outname = args['output']

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

# Setup tree directory to store output - make one if it doesn't exist

if not os.path.isdir('{0}/trees'.format(working_dir)):
    os.mkdir('{0}/trees'.format(working_dir))

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
    if not os.path.isfile('{0}/trees/{1}.nwk'.format(working_dir,
                          name)):
        tree = get_best_TC_tree(dv, gm, labels_file, tree_files, name)
        print tree
        tree.write_to_file('{0}/trees/{1}.nwk'.format(working_dir,
                           name), metadata=True)
    else:
        tree = Tree()
        tree.read_from_file('{0}/trees/{1}.nwk'.format(working_dir,
                            name))

    dv_matrix_strip_header = '\n'.join(dv_matrix.split('\n'
            )[2:]).rstrip()
    labels_strip_header = labels.split('\n')[1].rstrip()
    record = TCSeqRec()
    record.dv = [(dv_matrix_strip_header, labels_strip_header)]
    record.tree = tree
    record.name = name
    record.headers = labels_strip_header.split()
    record.sequences = ['' for _ in record.headers]
    record._update()
    records.append(record)

collection = SequenceCollection(records=records, get_distances=False,
                                gtp_path=os.environ['GTP_PATH'])
collection.put_distance_matrices('rf')
T = \
    collection.Clustering.run_spectral_rotate(collection.distance_matrices['rf'
        ])
collection.partitions[T] = Partition(T)
collection.clusters_to_partitions[('rf', 'spectral_rotate', max(T))] = T
collection.concatenate_records()
cluster_recs = collection.get_cluster_records()

number_of_clusters = len(cluster_recs)
for j in range(number_of_clusters):
    record = cluster_recs[j]
    record_dv = record.dv[0]
    labels = record.dv[1]

    # Write some temp files from our concatenated record
    # as input for tree collection -
    # ..._dv.txt     = concatenated distance matrices
    # ..._map.txt    = updated genome map - may have gained new
    #                  species in the concatenation, also labels
    #                  can change order
    # ..._labels.txt = updated labels list - in case concatenation
    #                  changes the order
    # NB: we can keep the initial guide trees from before

    tmpdir = '/nfs/nobackup/goldmans/kg'
    filename = record._write_temp_tc(tmpdir=tmpdir,
            make_guide_tree=False, hashname=True)
    tmp_dv = '{0}/{1}_dv.txt'.format(tmpdir, filename)
    tmp_map = '{0}/{1}_map.txt'.format(tmpdir, filename)
    tmp_labels = '{0}/{1}_labels.txt'.format(tmpdir, filename)

    tree = get_best_TC_tree(tmp_dv, tmp_map, tmp_labels, tree_files,
                            record.name)
    for tmpfile in [tmp_dv, tmp_map, tmp_labels]:
        os.remove(tmpfile)

    if not os.path.isfile('{0}/trees/cluster{1}.nwk'.format(working_dir,
                          j)):
        tree.write_to_file('{0}/trees/cluster{1}.nwk'.format(working_dir,
                           j), metadata=True)
    else:
        tree = Tree()
        tree.read_from_file('{0}/trees/cluster{1}.nwk'.format(working_dir,
                            j))
    record.tree = tree

    print j + 1
    print multiwordReplace(tree.newick, names)

result = cPickle.dump(collection,
                      open('{0}/{1}.pickle'.format(working_dir,
                      outname), 'w'))
