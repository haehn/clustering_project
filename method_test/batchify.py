#!/usr/bin/env python

import glob
import argparse
import re
import os
import math
import sys
from sequence_record import TCSeqRec

"""
Usage: python batchify.py -b <-path-/basedir-> -p <dna_alignments>
"""

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
parser.add_argument('-b', '--basedir', help='Path of directory containing all the simulation replicates', type=fpath, default='.')
parser.add_argument('-p', '--phylip_dir', help='Subpath of simdir in which to find the phylip files', type=fpath)
parser.add_argument('-r', '--just_read', help='Quit program after reading files, without writing to disk', action='store_true')
args = vars(parser.parse_args())
basedir = args['basedir']
nsims = len(glob.glob('{0}/sim*'.format(basedir)))
just_read = args['just_read']
#phydirs = [basedir + '/sim{0}/'.format(i) + args['phylip_dir'] for i in range(1,nsims+1)]
phydirs = [basedir + '/level{0}'.format(i) + '/sim{0}/'.format(j) + args['phylip_dir'] for i in range(1,11) for j in range(1,51)]
l=[]

print 'Reading locations of all phylip files...'

for d in phydirs:
    for f in glob.glob('{0}/*.phy'.format(d)):
        if not os.path.isfile(f[:f.rindex('.')]+'.ml.pickle'):
            print f
            l.append(os.path.abspath(f))
print '(done).'
if just_read: sys.exit()
print 'Making necessary batch directories...'
for i in range(int(math.floor((len(l)-1)/1000)+1)):
    dname = '{0}/{1}-{2}'.format(basedir,i*1000+1,(i+1)*1000)
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
    d = which_dir(i,basedir)
    with open('{0}/phylip.{1}'.format(d,i),'w') as file:
        file.write(p)
print 'Finished.'
