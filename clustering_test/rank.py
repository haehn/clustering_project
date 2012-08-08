#!/usr/bin/env python

import glob
import argparse
import os
from collections import defaultdict


def fpath(s):
    """
    Helper function used when passing filepath arguments with argparse module.
        Trims all '/' characters from the end of the path string.
    """

    while s.endswith('/'):
        s = s[:-1]
    return s


def quick_parse(d):
    getval = lambda x: x.rstrip().split()[-1]
    euc = glob.glob('{0}/euc*'.format(d))
    geo = glob.glob('{0}/geo*'.format(d))
    sym = glob.glob('{0}/sym*'.format(d))

    eucdic = {}
    geodic = {}
    symdic = {}

    for f in euc:
        with open(f) as file:
            for line in file:
                if line.startswith('Method'):
                    method = getval(line)
                if line.startswith('Score'):
                    score = float(getval(line))
                if line.startswith('Varinf'):
                    vi = float(getval(line))
            eucdic[method] = (score, vi)

    for f in geo:
        with open(f) as file:
            for line in file:
                if line.startswith('Method'):
                    method = getval(line)
                if line.startswith('Score'):
                    score = float(getval(line))
                if line.startswith('Varinf'):
                    vi = float(getval(line))
            geodic[method] = (score, vi)

    for f in sym:
        with open(f) as file:
            for line in file:
                if line.startswith('Method'):
                    method = getval(line)
                if line.startswith('Score'):
                    score = float(getval(line))
                if line.startswith('Varinf'):
                    vi = float(getval(line))
            symdic[method] = (score, vi)

    return (eucdic, geodic, symdic)


def get_rankings_by_score(dic):
    if not dic: return []
    get_score = lambda x, y: x[y][0]
    second_item = lambda z: z[1]
    scores = sorted([(key, get_score(dic, key)) for key in dic],
                    key=second_item, reverse=True)
    best = scores[0][1]
    eps = 0.01
    scorers = [x[0] for x in scores if abs(x[1] - best) < eps]
    num_scorers = len(scorers)
    return [(x, 1.0 / num_scorers) for x in scorers]


def get_rankings_by_vi(dic):
    if not dic: return []
    get_vi = lambda x, k: x[k][1]
    second_item = lambda z: z[1]
    vis = sorted([(key, get_vi(dic, key)) for key in dic],
                 key=second_item)
    best = vis[0][1]
    eps = 0.01
    scorers = [x[0] for x in vis if abs(x[1] - best) < eps]
    num_scorers = len(scorers)
    return [(x, 1.0 / num_scorers) for x in scorers]


parser = argparse.ArgumentParser(prog='rank.py')
parser.add_argument('-d', '--base-directory', type=fpath)
args = vars(parser.parse_args())

score_tracker = {'euc': defaultdict(float), 'geo': defaultdict(float),
                 'sym': defaultdict(float)}
vi_tracker = {'euc': defaultdict(float), 'geo': defaultdict(float),
              'sym': defaultdict(float)}

basedir = args['base_directory']
simdirs = [x for x in glob.glob('{0}/sim*'.format(basedir)) if os.path.isdir(x)]

for directory in simdirs:
    print directory
    (e, g, s) = quick_parse(directory)
    for (method, rankscore) in get_rankings_by_score(e):
        score_tracker['euc'][method] += rankscore
    for (method, rankscore) in get_rankings_by_score(g):
        score_tracker['geo'][method] += rankscore
    for (method, rankscore) in get_rankings_by_score(s):
        score_tracker['sym'][method] += rankscore
    for (method, rankscore) in get_rankings_by_vi(e):
        vi_tracker['euc'][method] += rankscore
    for (method, rankscore) in get_rankings_by_vi(g):
        vi_tracker['geo'][method] += rankscore
    for (method, rankscore) in get_rankings_by_vi(s):
        vi_tracker['sym'][method] += rankscore

print 'Results by score:'
print
print 'Euclidean distance'
for k in sorted(score_tracker['euc']):
    print '{0} {1}'.format(k,score_tracker['euc'][k])
print
print 'Geodesic distance'
for k in sorted(score_tracker['geo']):
    print '{0} {1}'.format(k,score_tracker['geo'][k])
print
print 'Robinson-Foulds distance'
for k in sorted(score_tracker['sym']):
    print '{0} {1}'.format(k,score_tracker['sym'][k])

print
print 'Results by Variation of Information:'
print
print 'Euclidean distance'
for k in sorted(vi_tracker['euc']):
    print '{0} {1}'.format(k,vi_tracker['euc'][k])
print
print 'Geodesic distance'
for k in sorted(vi_tracker['geo']):
    print '{0} {1}'.format(k,vi_tracker['geo'][k])
print
print 'Robinson-Foulds distance'
for k in sorted(vi_tracker['sym']):
    print '{0} {1}'.format(k,vi_tracker['sym'][k])
