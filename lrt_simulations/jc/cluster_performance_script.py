#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import cPickle
import os
import sys
from sequence_collection import SequenceCollection
from clustering import Clustering
from tree import Tree
from sequence_record import TCSeqRec
from result import Result
import glob
import dendropy as dpy
import numpy as np
import re
import shutil


# Command-line arguments handled by argparse module
parser = argparse.ArgumentParser(prog='clustersuccess_script.py',
                description='Check success of clustering procedure'
                )
parser.add_argument(
    '-b',
    '--basefile',
    dest='base',
    help='path to base treeset file',
    type=str,
    default=''
    )
parser.add_argument(
    '-o',
    '--outfile',
    dest='out',
    help='output pickle file',
    type=str,
    default=''
    )
parser.add_argument(
    '-s',
    '--set',
    dest='set',
    help='choice of treeset',
    type=int,
    default='1'
    )

args = vars(parser.parse_args())
choice = args['set']
### END ARGPARSE BIT

# Essential files
base_file = args['base']
helper = os.environ['DARWINHELPER']
output_file = args['out']
tmpdir = os.environ['TEMPORARY_DIRECTORY']

for f in [base_file, helper]:
    if not os.path.isfile(f):
        print 'Had a problem finding this file: {0}'.format(f)
        sys.exit(0)
### END FILE CHECKS

# main

trees = open(base_file).readlines()
treesets = []
while trees:
    treeset = []
    while not trees[0] == '\n':
        treeset.append(trees[0])
        trees = trees[1:]
    treesets.append(treeset)
    trees=trees[1:]

treesets = treesets[choice-1:choice] ### REMOVE THIS LINE
truth = np.array( [1]*10+[2]*10+[3]*10 )

for i,treeset in enumerate(treesets):
    treeset = treesets[i]
    dm = Clustering().get_distance_matrix( trees = [Tree(newick=tree).pam2sps('pam2sps') for tree in treeset], metric='sym' )
    avdist_sym = np.mean(dm[np.triu_indices(len(dm),1)])
    dm = Clustering().get_distance_matrix( trees = [Tree(newick=tree).pam2sps('pam2sps') for tree in treeset], metric='euc' )
    avdist_euc = np.mean(dm[np.triu_indices(len(dm),1)])
    dm = Clustering().get_distance_matrix( trees = [Tree(newick=tree).pam2sps('pam2sps') for tree in treeset], metric='rf' )
    avdist_rf = np.mean(dm[np.triu_indices(len(dm),1)])
    dm = Clustering().get_distance_matrix( trees = [Tree(newick=tree).pam2sps('pam2sps') for tree in treeset], metric='geodesic' )
    avdist_geo = np.mean(dm[np.triu_indices(len(dm),1)])
    msa_dir = '{0}/exp{1}_{2}'.format(tmpdir,i,choice)
    if not os.path.isdir(msa_dir):
        os.mkdir(msa_dir)
    number = 0
    for j, tree in enumerate(treeset):
        tmptree = '{0}/tmptree{1}_{2}_{3}.nwk'.format(tmpdir,i,j,choice)
        params = '{0}/{1}_{2}_{3}_params.drw'.format(tmpdir,i,j,choice)
        with open(tmptree,'w') as file:
            file.write(tree)
        Result().write_ALF_parameters(
            'alfsim{0}_{1}_{2}'.format(i,j,choice),
            tmpdir,
            'alftmp{0}_{1}_{2}'.format(i,j,choice),
            10,
            250,
            tmptree,
            params
            )
        os.system('alfsim {0}'.format(params))
        os.remove(params)
        os.remove(tmptree)
        alf_newick = \
            open('{0}/alftmp{1}_{2}_{3}/alfsim{1}_{2}_{3}/RealTree.nwk'.format(tmpdir,
                 i,j,choice)).read()
        replacement_dict = dict(zip(re.findall(r'(\w+)(?=:)',
                                alf_newick),
                                re.findall(r'(\w+)(?=:)',
                                tree)))
        print alf_newick
        print tree
        print replacement_dict
        fasta_files = glob.glob('{0}/alftmp{1}_{2}_{3}/alfsim{1}_{2}_{3}/MSA/*dna.fa'.format(tmpdir,
                 i,j,choice))
        for fasta_file in fasta_files:
            number += 1
            name = 'gene{0:0>3}'.format(number)
            record = TCSeqRec(fasta_file)
            new_headers = [replacement_dict[x[:x.rindex('/')]] for x in
                       record.headers]
            print new_headers
            sequences = record.sequences
            new_record = TCSeqRec(headers=new_headers, sequences=sequences,
                                    name=name)
            new_record.write_fasta(outfile='{0}/{1}.fas'.format(msa_dir,name))
        shutil.rmtree('{0}/alftmp{1}_{2}_{3}'.format(tmpdir,
                 i,j,choice))
    sim = SequenceCollection(msa_dir, datatype='dna',helper=os.environ['DARWINHELPER'])
    sim.put_trees_parallel(program='phyml',model='JC69',datatype='nt',ncat=1)
    sim.put_partitions(metrics=['sym','rf','euc','geodesic'], linkages=['single','complete','average','ward'], nclasses=[3])
    sim.put_clusters()
    sim.put_cluster_trees_parallel(program='phyml',model='JC69',datatype='nt',ncat=1)
    d = sim.get_clusters()
    print 'Av. Dist\t(sym):\t{0}'.format(avdist_sym)
    print 'Av. Dist\t(euc):\t{0}'.format(avdist_euc)
    print 'Av. Dist\t(rf):\t{0}'.format(avdist_rf)
    print 'Av. Dist\t(geo):\t{0}'.format(avdist_geo)
    with open(output_file,'w') as file:
        for k in sorted(d):
            partition = sim.get_partitions()[k]
            vi = sim.clustering.variation_of_information(partition,truth)
            score = d[k].score
            if k[0] == 'sym':
                file.write('sym\t{0}\t{1}\t{2}\t{3}'.format(k[1],avdist_sym,vi,score))
            elif k[0] == 'rf':
                file.write('rf\t{0}\t{1}\t{2}\t{3}'.format(k[1],avdist_rf,vi,score))
            elif k[0] == 'geodesic':
                file.write('geodesic\t{0}\t{1}\t{2}\t{3}'.format(k[1],avdist_geo,vi,score))
            elif k[0] == 'euc':
                file.write('euc\t{0}\t{1}\t{2}\t{3}'.format(k[1],avdist_euc,vi,score))
