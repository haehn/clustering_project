#!/usr/bin/env python

import glob
import os
import re
import sys

indir = sys.argv[1]
dirs = glob.glob('{0}/level*'.format(indir))
for d in dirs:
    if not os.path.isdir(d):
        dirs.remove(d)

r=re.compile(r'(?<=level)[\d+]+')
outf = sys.argv[2]
P = sys.argv[3]
with open(outf, 'a') as f:
    for d in dirs:
        print d
        try:
            level = r.search(d).group()
        except AttributeError:
            sys.exit(1)
        with open('{0}/summary_scores.txt'.format(d)) as data:
            for line in data:
                line=line.rstrip()
                dist, meth, lnl, delta, vi = line.split(',')
                f.write(','.join([str(x) for x in [P, level, dist, meth, lnl, delta, vi]])+'\n')

