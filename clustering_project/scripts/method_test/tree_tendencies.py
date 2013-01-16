#!/usr/bin/env python

import glob
import re

base = '/nfs/research/goldman/kevin/regime2'
names_dict = dict(weighted_spr='SPR', newcoal='Coalescent', nni='NNI')

topper = lambda x: x[x.rindex('/') + 1:]
number_grabber = lambda x: re.search(r'(\d+)', x).group()


def extractor(f):
    lspl = lambda l: l.rstrip().split('\t')[:2]
    rec = False
    d = {}
    with open(f) as inf:
        for line in inf:
            if line.startswith('Class base tree distances'):
                rec = True
            elif line.startswith('Geodesic class stats'):
                rec = False
                break
            if (line.startswith('geodesic\t')
                or line.startswith('euclidean\t')
                or line.startswith('RF\t') or line.startswith('wRF\t')) \
                and rec:
                (k, v) = lspl(line)
                d[k] = v
    return d


def write_line(
    gen='',
    lev='',
    rep='',
    met='',
    dst='',
    hea=False,
    ):
    d = dict(geodesic='Geodesic', euclidean='Euclidean',
             RF='Robinson-Foulds', wRF='Weighted Robinson-Foulds')
    if hea:
        return 'GENERATOR\tLEVEL\tDISTANCE\tMETHOD\tREP\n'
    return '{0}\t{1}\t{2}\t{3}\t{4}\n'.format(gen, lev, dst, d[met], rep)

with open('tree_tendencies.tsv','w') as outf:
    outf.write(write_line(hea=True))
    for d in glob.glob(base + '/*'):
        if topper(d) in names_dict:
            print d
            generator = names_dict[topper(d)]
            for leveldir in glob.glob('{0}/level*'.format(d)):
                level = number_grabber(topper(leveldir))
                for simdir in glob.glob('{0}/sim*'.format(leveldir)):
                    rep = number_grabber(topper(simdir))
                    dists = \
                        extractor('{0}/treedistances.txt'.format(simdir))
                    for (m,d) in dists.items():
                        outf.write(write_line(generator, level, rep, m, d))

