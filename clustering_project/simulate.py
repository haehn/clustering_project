#!/usr/bin/env python

from seqsim import SeqSim
import numpy as np
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


help_regime = \
    '''
The data can be simulated according to several regimes:

1 - All genes come from the same tree
-------------------------------------

The genes will all be evolved along an identical tree.
This is the same as setting the number of classes to 1.


2 - Genes come from more than one tree
--------------------------------------

The genes are distributed into K classes. The genes in each class
are evolved along their own tree. The similarity between class
trees can be controlled.


3 - Genes come from more than one tree, and evolutionary rates
    vary between genes
--------------------------------------------------------------

The genes are distributed into K classes. The genes in each class
are evolved along the same tree topology, but the branch lengths
are scaled by a parameter for each gene.


4 - Genes come from more than one tree, and evolutionary rates
    vary between genes and lineages
----------------------------------------

The genes are distributed into K classes. The genes in each class
are evolved along the same tree topology, but branch lengths
are free to vary within each tree for each gene.

'''

help_tuning = \
    '''
Tuning parameter for adjusting assignment of genes into
classes. High values produce skewed assignments. 0 produces
a balanced assignment.

'''

help_permutations = \
    '''
The simulator generates a master tree according
to some process - random branching, Yule Pure Birth,
coalescence. The class trees are derived from the master
tree by applying some kind of permutation to it. 

This parameter controls the extent of the permutation
applied to the master tree when producing class trees.

It is the integer number of operations to apply when using 
nearest-neighbour interchanges or subtree prune and regrafts, 
or the length in coalescent units to scale the tree to when 
generating gene trees from a species tree - longer trees 
produce more congruent gene trees

'''

help_master = \
    '''
Method for producing the master tree.
Valid options are:

1 - random_topology
-------------------

The tree topology is produced by producing a node and 
randomly connecting it to an existing edge.
Branch lengths are drawn randomly from a gamma distribution

2 - random_yule
---------------

An ultrametric speciation tree is produced using a uniform
pure-birth process.

3 - random_coal
---------------

An ultrametric coalescent tree is produced by randomly 
coalescing leaf nodes.

'''

help_permuter = \
    '''
Operation that produces class trees from the master tree

1 - nni
-------
Use nearest-neighbour interchanges to create class trees

2 - spr
-------
Use subtree-prune-and-regraft moves to create class trees.
This should be used on ultrametric trees to represent
lateral gene transfer

3 - coal
--------
Sample a gene tree from the master (species) tree using
the constrained Kingman approach

'''

help_unique = \
    '''
If selected, constrains class trees to have different topologies

'''
desc = 'Run ALF simulator to generate different topological classes'
parser = argparse.ArgumentParser(prog='testsim.py', description=desc,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-k', '--classes', help='Number of classes',
                    type=int, default=2)
parser.add_argument('-n', '--species', help='Number of species',
                    type=int, default=20)
parser.add_argument(
    '-m',
    '--genes',
    help='Number of genes',
    type=int,
    nargs='+',
    default=[25, 25],
    )
parser.add_argument('-r', '--regime', help=help_regime, type=int,
                    default=2)
parser.add_argument('-t', '--tune', help=help_tuning, type=int,
                    default=None)
parser.add_argument('-master', help=help_master, type=str,
                    default='random_yule')
parser.add_argument('-c', '--class_permuter', help=help_permuter,
                    type=str, choices=['nni','spr', 'coal'], default='nni')
parser.add_argument('-p', '--permutations', help=help_permutations,
                    type=int, default=1)
parser.add_argument('-d', '--directory', help='Base output directory\n'
                    , type=fpath, default='.')
parser.add_argument('-g', '--geodesic',
                    help='path to gtp.jar, used to calculate geodesic distance\n'
                    , type=fpath, default='./class_files')
parser.add_argument('-tmp', '--temp-directory',
                    help='Directory to use for temp files\n',
                    type=fpath, default='/tmp')
parser.add_argument('-i', '--indels',
                    help='Simulate indels (default=no)\n',
                    action='store_true')
parser.add_argument('-ratevar', '--ratevar',
                    help='Simulate rate variation among sites (default=no)\n'
                    , action='store_true')
parser.add_argument('-q', '--quiet', help='Less printing to screen\n',
                    action='store_true')
parser.add_argument('-u', '--unique', help=help_unique,
                    action='store_true')
args = vars(parser.parse_args())

K = args['classes']
mk = args['genes']
if len(mk) > 1:
    M = sum(mk)
    K = len(mk)
else:
    M = mk[0]
    K = 1
    mk = None

n = args['species']
tune = args['tune']
regime = args['regime']
permutation_strength = args['permutations']
permutation_type = args['class_permuter']
master_tree_generator = args['master']
filepath = args['directory']
if 'TEMPORARY_DIRECTORY' in os.environ:
    tmpdir = os.environ['TEMPORARY_DIRECTORY']
else:
    tmpdir = args['temp_directory']
indels = args['indels']
ratevar = args['ratevar']
quiet = args['quiet']
gtp_path = args['geodesic']

sim = SeqSim('base')

sim.hky_model(
    alpha=3.551,
    beta=1,
    Afreq=0.31,
    Cfreq=0.18,
    Gfreq=0.21,
    Tfreq=0.3,
    )

if indels:
    sim.indels()

if ratevar:
    sim.rate_variation()

print 'Regime = {0}'.format(regime)
print 'Parameter settings:'
params = '\n  '.join([
    '  K: {0}'.format(K),
    'M: {0}'.format(M),
    'mk: {0}'.format(mk),
    'n: {0}'.format(n),
    'tune: {0}'.format(tune),
    'permutation type: {0}'.format(permutation_type),
    'permutation strength: {0}'.format(permutation_strength),
    'master tree generator: {0}'.format(master_tree_generator),
    'filepath: {0}'.format(filepath),
    'tmpdir: {0}'.format(tmpdir),
    'gtp path: {0}'.format(gtp_path),
    'indels: {0}'.format(indels),
    'ratevar: {0}'.format(ratevar),
    ])
print params
print '''


'''
sim.simulate_set(
    K=K,
    M=M,
    n=n,
    tune=tune,
    regime=regime,
    branch_length_func=np.random.gamma,
    inner_edge_params=(3.2, 0.029),
    leaf_params=(2.2, 0.097),
    scale_func=np.random.gamma,
    mk=mk,
    master_tree_generator_method=master_tree_generator,
    class_tree_permuter=permutation_type,
    guarantee_unique=True,
    num_permutations=permutation_strength,
    scale_params=(2, 0.5),
    gene_length_kappa=6,
    gene_length_theta=83.33,
    gene_length_min=50,
    filepath=filepath,
    tmpdir=tmpdir,
    gtp_path=gtp_path,
    unit_is_pam=True,
    quiet=quiet,
    )
print '(Done).'
