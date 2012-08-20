#!/usr/bin/env python

from sequence_record import TCSeqRec
import argparse
import os
import cPickle
import sys

# Runs a batch of phyml jobs using LSB JobArrays. Batches are prepared w/
# batchify.py
# Usage: bsub -J "jobname[#-##]" python runbatchphyml.py -f <path>/phylip.

# bash shell expansion completes the input filename with the jobindex which is
# determined in the jobname[#-##] flag - so make sure these match

def fpath(s):
    """
    Helper function used when passing filepath arguments with argparse module.
        Trims all '/' characters from the end of the path string.
    """
    while s.endswith('/'):
        s = s[:-1]
    return s

getname = lambda x: x[x.rindex('/')+1:x.rindex('.')]


parser = argparse.ArgumentParser(prog='phymlwrap.py')
parser.add_argument('-f', '--infile', help='input file', type=fpath, default='.')
args = vars(parser.parse_args())
try: 
    index = os.environ['LSB_JOBINDEX']
except:
    index = None
infile = args['infile']
if index:
    infile += index

if not infile[-1].isdigit():
    print '{0} is not correct'.format(infile)
    sys.exit(1)

if not os.path.isfile(infile):
    print 'Input file not found:\n{0}'.format(os.path.abspath(infile))
    sys.exit(0)

with open(infile) as file:
    target = os.path.abspath(file.read())

if not os.path.isfile(target):
    print 'Target file not found:\n{0}'.format(os.path.abspath(target))
    sys.exit(2)

print target
parent_dir = os.path.dirname(target)
name = getname(target)

#if os.path.isfile('{0}/{1}.nj.pickle'.format(parent_dir, name)):
#    os.remove('{0}/{1}.nj.pickle'.format(parent_dir, name))
#if os.path.isfile('{0}/{1}.pickle'.format(parent_dir, name)):
#    os.remove('{0}/{1}.pickle'.format(parent_dir, name))
record = TCSeqRec(target, file_format='phylip', name=name)
if 'phyml_clustering' in target:
    record.get_phyml_tree(model='GTR',ncat=4,datatype='nt')
elif 'bionj_clustering' in target:
    record.get_bionj_tree(model='GTR',ncat=4,datatype='nt')
cPickle.dump(record, open('{0}/{1}.pickle'.format(parent_dir, name),'w'))
#record.get_bionj_tree(model='GTR',ncat=4,datatype='nt')
#cPickle.dump(record, open('{0}/{1}.nj.pickle'.format(parent_dir, name),'w'))
#record = TCSeqRec(target, file_format='phylip', name=name)
#record.get_phyml_tree(model='GTR',ncat=4,datatype='nt')
#cPickle.dump(record, open('{0}/{1}.pickle'.format(parent_dir, name), 'w'))
