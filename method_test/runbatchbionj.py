#!/usr/bin/env python

import argparse
from sequence_record import TCSeqRec
import cPickle
import os
import glob
import re
import sys

"""
Usage: bsub -o /dev/null -e error.%J.%I -J "jobname[1-x]" bash <path>/tempdir_wrapper.sh python doclustering.py -d <path>/nni/level1/sim_

Where the path endswith an underscore above, JobArray will fill in the number of the sim from 1 to x

"""

def fpath(s):
    """
    Helper function used when passing filepath arguments with argparse module.
        Trims all '/' characters from the end of the path string.
    """
    while s.endswith('/'):
        s = s[:-1]
    return s

sort_key = lambda item: tuple((int(num) if num else alpha) for (num,alpha) in re.findall(r'(\d+)|(\D+)', item))
getname = lambda x: x[x.rindex('/')+1:x.rindex('.')]

parser = argparse.ArgumentParser(prog='doclustering.py')
parser.add_argument('-d', '--directory', help='input directory', type=fpath, default='.')
parser.add_argument('-p', '--phylip_dir', help='Subpath of simdir in which to find the phylip files', type=fpath)
args = vars(parser.parse_args())
index = os.environ['LSB_JOBINDEX']
indir = os.path.abspath(args['directory'])
if index and index != '0':
    indir += index
print 'Working on {0}'.format(indir)
phylip_dir = args['phylip_dir']
try:
    tmpdir = os.environ['TEMPORARY_DIRECTORY']
except:
    tmpdir = '/tmp'

working_dir = '/'.join([indir, phylip_dir])
phylip_files = sorted(glob.glob('{0}/*.phy'.format(working_dir)), key=sort_key)
print 'Working on {0}'.format(working_dir)

for f in phylip_files:
    name = getname(f)
    print 'Getting BIONJ tree for {0}'.format(name)
    seqrec = TCSeqRec(f, file_format='phylip', name=name)
    seqrec.datatype='dna'
    seqrec.get_bionj_tree(model='GTR', ncat=4, datatype='nt', tmpdir=tmpdir)
    cPickle.dump(seqrec, open('{0}/{1}.nj.pickle'.format(working_dir, name),'w'))# In future let's just pickle trees; sequences already stored on disk
