#!/usr/bin/env python

import glob
import os
import re
import sys
import numpy as np
from collections import defaultdict

indir = sys.argv[1]
simdirs = glob.glob('{0}/sim*'.format(indir))
assert len(simdirs) == 50

dist_data = {
    'rf': defaultdict(list),
    'euc': defaultdict(list),
    'geo': defaultdict(list),
}
# method_data = {
#     'single': defaultdict(list),
#     'complete': defaultdict(list),
#     'average': defaultdict(list),
#     'ward': defaultdict(list),
#     'kmedoids': defaultdict(list),
#     'spectral': defaultdict(list),
#     'MDS': defaultdict(list),
# }
regex = re.compile('(?<=level)[\d+]+')
for d in simdirs:
    try:
        with open('{0}/scores.txt'.format(d)) as f:
            
            lines = f.read().split('\n')[:-1]
            trueLnL = float(lines.pop().split(',')[1])
            for line in lines:
                line = line.split(',')
                dist,method = line[0].split()
                lnl = float(line[1])
                adj_lnl = lnl - trueLnL
                vi = float(line[2])
                dist_data[dist][method].append((lnl, adj_lnl, vi))
                #method_data[method][dist].append((lnl, adj_lnl, vi))
    except IOError:
        #i = regex.search(d).group(0)
        print d
        #os.system( 'rm -r {0}'.format(d) )
        #os.system('bsub -o /dev/null ./tempdir_wrapper.sh python simulate.py  -k 4 -n 20 -m 60 -r 2 -t 0 -master random_yule -c coal -p {0} -d {1} -g class_files -tmp=/tmp/kg001 -i -ratevar -u'.format(i, d))
        #os.system('bsub -o /dev/null ../tempdir_wrapper.sh python calc_scores.py -d {0} -s dna -n 4'.format(d))
        #print 'Can\'t open file {0}/scores.txt'.format(d)
        # sys.exit(1)
        
with open('{0}/summary_scores.txt'.format(indir), 'w') as outf:
    for k in dist_data:
        for m in dist_data[k]:
            data = dist_data[k][m]
            lnl_bar = np.mean([x[0] for x in data])
            delta_bar = np.mean([x[1] for x in data])
            VI_bar = np.mean([x[2] for x in data])
            outf.write(','.join( [str(x) for x in [k, m, lnl_bar, delta_bar, VI_bar]] )+'\n')

