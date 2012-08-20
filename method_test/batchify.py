#!/usr/bin/env python

import glob
import argparse
import re
import os
import math
import sys
from sequence_record import TCSeqRec

def which_dir(n,base):
    avail_dirs = [x for x in glob.glob('{0}/*'.format(base)) if os.path.isdir(x) and x[x.rindex('/')+1].isdigit()]
    for x in avail_dirs:
        lower,upper = x[x.rindex('/')+1:].split('-')
        if int(lower) <= n <= int(upper):
            return x

def fpath(s):
    """
    Helper function used when passing filepath arguments with argparse module.
        Trims all '/' characters from the end of the path string.
    """
    while s.endswith('/'):
        s = s[:-1]
    return s

sort_key = lambda item: tuple((int(num) if num else alpha) for (num,alpha) in re.findall(r'(\d+)|(\D+)', item))

parser = argparse.ArgumentParser(prog='batchify.py')
parser.add_argument('-d', '--directory', help='input directory', type=fpath, default='.')
args = vars(parser.parse_args())
regdir = args['directory']
#basedirs = glob.glob('{0}/nni*'.format(regdir))
basedirs = ['{0}/nni{1}'.format(regdir,i) for i in range(2,6)]
l=[]

print 'Reading locations of all phylip files...'
for d in basedirs:
    for sd in [x for x in glob.glob('{0}/*'.format(d)) if os.path.isdir(x)]:
        for f in glob.glob('{0}/phyml_clustering/*.phy'.format(sd))+glob.glob('{0}/bionj_clustering/*.phy'.format(sd)):
            if not os.path.isfile(f[:-3]+'pickle'):# or not os.path.isfile(f[:-3]+'nj.pickle'):
                #print f
                l.append(os.path.abspath(f))
print '(done).'
print l
print 'Making necessary batch directories...'
for i in range(int(math.floor(len(l)/1000)+1)):
    dname = '{0}/{1}-{2}'.format(regdir,i*1000+1,(i+1)*1000)
    if not os.path.isdir(dname):
        print '   ',dname
        os.mkdir(dname)
print '(done).'

if not all(l):
    print 'None values found in location list, quitting.'
    sys.exit(0)

print 'Writing locations to batch directories...'
for i,p in enumerate(l,start=1):
    if i%1000 == 0:
        print i
    d = which_dir(i,regdir)
    with open('{0}/phylip.{1}'.format(d,i),'w') as file:
        file.write(p)
print 'Finished.'
