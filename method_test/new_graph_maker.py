#!/usr/bin/env python
print 'importing...'
import argparse
import glob
import numpy as np
from pylab import *
from collections import defaultdict
import re
import sys
print 'done.'

def fpath(s):
    while s.endswith('/'):
        s = s[:-1]
    return s

print 'parsing commandline...'
parser = argparse.ArgumentParser(prog='graph_maker.py')
parser.add_argument('-d', '--directory', help='Input directory', type=fpath, default='.')
parser.add_argument('-m', '--metric', help='Distance metric', type=str, default='sym')
parser.add_argument('-s', '--score', help='Either <\'lnl\'> or <\'vi\'>', type=str, default='lnl')
parser.add_argument('-p', '--program', help='Either <\'phyml\'> or <\'bionj\'>', type=str, default='phyml')
parser.add_argument('-x', '--xaxis', help='Which xaxis label to use: nni, spr or coal', type=str, default='nni')
parser.add_argument('-o', '--outfile', help='Outfile (for pdf image)', type=fpath, default='plot.pdf')
parser.add_argument('-n', '--normalise', action='store_true')
args = vars(parser.parse_args())
indir = args['directory']
metric = args['metric']
scoring_system = args['score']
program = args['program']
outfile = args['outfile']
xax = args['xaxis']
normalise = args['normalise']
subdirs = glob.glob('{0}/level*'.format(indir))
print 'done.'
getlevel = lambda x: re.search(r'[\d+]+', x[x.rindex('/')+1:]).group()

scoresdic = defaultdict(dict)

if not program in ['phyml', 'bionj']:
    print 'Unrecognised program selected:', program
    sys.exit(-1)

for d in subdirs:
    level = int(getlevel(d))
    with open('{0}/{1}_{2}_results.txt'.format(d, program, scoring_system)) as infile:
        for line in infile:
            name, l = line.rstrip().split('\t')
            mean_score = np.mean([float(x) for x in l.split(' ')])
            scoresdic[level][name] = [float(x) for x in l.split(' ')]

if scoring_system == 'lnl':
    ylabel_text = 'Difference in log-likelihood from true partitioning'
elif scoring_system == 'vi':
    ylabel_text = 'Partitioning distance (Variation of Information)'

if metric == 'euc':
    title_text = 'Euclidean distances'
elif metric == 'sym':
    title_text = 'Robinson-Foulds distances'
elif metric == 'geo':
    title_text = 'Geodesic distances'
else: title_text = 'All distances'

if program == 'phyml':
    suptitle_text = 'Maximum Likelihood trees'
elif program == 'bionj':
    suptitle_text = 'BioNJ trees'

print 'setting up chart...'
fig = figure()
ax = fig.add_subplot(111)
ax.set_xlim( [0, max(scoresdic.keys())+1])
if xax == 'nni':
    xlabel('Number of nnis between classes')
elif xax == 'spr':
    xlabel('Number of sprs between classes')
elif xax == 'coal':
    xlabel('Species tree depth (coalescent units)')
else:
    xlabel('')
ylabel(ylabel_text)
suptitle(suptitle_text)
title(title_text)

colors = dict(MDS = 'k',
            single = 'r',
            complete = 'b',
            average = 'y',
            ward = 'g',
            spectral = 'm',
            kmedoids = 'c')

markers = dict(euc = 'x',
    sym = (6,2,1),
    geo = '+')

linestyles = dict(euc = '--',
    sym = '-.',
    geo = ':')
print 'plotting points...'


xarray = sorted(scoresdic.keys())
yarrays = defaultdict(list)

for xval in xarray:    
    for key in scoresdic[xval]:
        score = scoresdic[xval][key]
        truescore = scoresdic[xval]['true']
        if normalise:
            yarrays[key].append(((np.mean(score)-np.mean(truescore))/np.mean(truescore)))
        else:
            yarrays[key].append((np.mean(score)-np.mean(truescore)))


for k in yarrays:
    if k == 'true': continue
    elif not metric in k and metric != 'all':
        continue
    yarray = yarrays[k]
    markerkey, colorkey = k.split(' ') 
    ax.plot(xarray,yarray, ls=linestyles[markerkey], lw=1, color=colors[colorkey], zorder=1)
    ax.scatter(xarray,yarray, marker=markers[markerkey], color=colors[colorkey], zorder=2)



# leg = ax.legend(loc='best')
fig.savefig(outfile)
