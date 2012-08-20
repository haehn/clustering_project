#!/usr/bin/env python

import argparse
import glob
import numpy as np
from pylab import *

def fpath(s):
    while s.endswith('/'):
        s = s[:-1]
    return s

parser = argparse.ArgumentParser(prog='graph_maker.py')
parser.add_argument('-d', '--directory', help='Input directory', type=fpath, default='.')
parser.add_argument('-m', '--metric', help='Distance metric', type=str, default='sym')
parser.add_argument('-s', '--score', help='Either <\'likelihood\'> or <\'vi\'>', type=str, default='likelihood')
parser.add_argument('-p', '--program', help='Either <\'phyml\'> or <\'bionj\'>', type=str, default='phyml')
parser.add_argument('-o', '--outfile', help='Outfile (for pdf image)', type=fpath, default='plot.pdf')
args = vars(parser.parse_args())
indir = args['directory']
metric = args['metric']
scoring_system = args['score']
program = args['program']
outfile = args['outfile']
nnidirs = glob.glob('{0}/nni*'.format(indir))

scoresdic = {}

for nnidir in nnidirs:

    name = nnidir[nnidir.rindex('/')+1:]
    number = int(name.split('nni')[1])
    
    if program == 'phyml':
        results_file = '{0}/phymlresults.txt'.format(nnidir)  
    elif program == 'bionj':
        results_file = '{0}/bionjresults.txt'.format(nnidir)
        
    results = [x.split('\t') for x in \
                open(results_file).read().split('\n')[:-1]]
    
    for i in range(len(results)):
        for j in range(1,len(results[i])):
            results[i][j] = float(results[i][j])
    
    scoresdic[number]=results

if metric == 'euc':
    title_text = 'Euclidean distances'
    rows = [1,2,3,4,5,6,7]
elif metric == 'sym':
    title_text = 'Robinson-Foulds distances'
    rows = [15,16,17,18,19,20,21]
elif metric == 'geo':
    title_text = 'Geodesic distances'
    rows = [8,9,10,11,12,13,14]

if program == 'phyml':
    suptitle_text = 'Maximum Likelihood trees'
elif program == 'bionj':
    suptitle_text = 'BioNJ trees'

if scoring_system == 'likelihood':
    column = 1
    ylabel_text = 'Difference in log-likelihood from true partitioning'
elif scoring_system == 'vi':
    column = 3
    ylabel_text = 'Partitioning distance (Variation of Information)'

x = sorted(scoresdic.keys())
methods = ['true'] + 3*['MDS','Average-linkage', 'Complete-linkage', 
                        'k-medoids', 'Single-linkage', 'Spectral',
                        'Ward\'s method']
fig = figure()
ax = fig.add_subplot(111)
xlabel('Number of NNIs between classes')
ylabel(ylabel_text)
suptitle(suptitle_text)
title(title_text)

for row in rows:
    y = []
    legend_text = methods[row]
    for key in x:
        true_value = scoresdic[key][0][column]
        observed_value = scoresdic[key][row][column]          
        y.append(observed_value - true_value)
    ax.plot(x,y,label=legend_text)

leg = ax.legend(loc='best')

fig.savefig(outfile)
    