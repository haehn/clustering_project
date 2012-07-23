#!/usr/bin/env python
# -*- coding: utf-8 -*-

from seqsim import SeqSim
import numpy as np
import argparse

def fpath(s):
    """
    Helper function used when passing filepath arguments with argparse module.
        Trims all '/' characters from the end of the path string.
    """
    while s.endswith('/'):
        s = s[:-1]
    return s

parser = argparse.ArgumentParser(prog='testsim.py',
                                 description='Run ALF simulator to generate different topological classes'
                                 )
parser.add_argument('-k', '--classes', help='Number of classes',
                    type=int, default=2)
parser.add_argument('-n', '--species', help='Number of species',
                    type=int, default=20)
parser.add_argument('-m', '--genes', help='Number of genes', type=int,
                    default=100)
parser.add_argument('-t', '--tune',
                    help='Tuning parameter for adjusting assignment of genes into classes'
                    , type=int, default=1)
parser.add_argument('-r', '--regime', help='Regime', type=int,
                    default=2)
parser.add_argument('-nni', '--nni', help='Number of NNIs', type=int,
                    default=0)
parser.add_argument(
    '-d',
    '-dir',
    '--directory',
    help='Base output directory',
    type=fpath,
    default='.',
    )
parser.add_argument('-tmp', '--temp-directory',
                    help='Directory to use for temp files', type=fpath,
                    default='/tmp')
parser.add_argument('-i', '--indels',
                    help='Simulate indels (default=no)',
                    action='store_true')
parser.add_argument('-ratevar', '--ratevar',
                    help='Simulate rate variation among sites (default=no)'
                    , action='store_true')
parser.add_argument('-q', '--quiet',
                    help='Less printing to screen'
                    , action='store_true')
args = vars(parser.parse_args())

K = args['classes']
M = args['genes']
n = args['species']
tune = args['tune']
regime = args['regime']
nni = args['nni']
filepath = args['directory']
tmpdir = args['temp_directory']
indels = args['indels']
ratevar = args['ratevar']
quiet = args['quiet']

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

sim.simulate_set(
    K=K,
    M=M,
    n=n,
    tune=tune,
    regime=regime,
    branch_length_func=np.random.gamma,
    inner_edge_params=(2.13, 0.029),
    leaf_params=(1.45, 0.097),
    gene_length_kappa=5.53,
    gene_length_theta=72.35,
    gene_length_min=10,
    filepath=filepath,
    scale_func=np.random.gamma,
    scale_params=(1, 1),
    nni=nni,
    tmpdir=tmpdir,
    unit_is_pam=True,
    quiet=quiet,
    )
print '(Done).'
