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

parser = argparse.ArgumentParser(prog='lsf_phyml.py')
parser.add_argument('-f', '--infile', help='input file', type=fpath, default='.')
parser.add_argument('-m', '--model', help='model', type=str, default='GTR')
parser.add_argument('-n', '--ncat', help='number of categories of gamma-distributed rate variation', type=int, default=4)
parser.add_argument('-d', '--datatype', help='datatype = nt (nucleotide), or aa (amino acid)', type=str, default='nt')
args = vars(parser.parse_args())

infile = args['infile'] + os.environ['LSB_JOBINDEX']
model = args['model']
datatype = args['datatype']
ncat = args['ncat']

if not os.path.isfile(infile):
    print 'Input file not found:\n{0}'.format(os.path.abspath(infile))
    sys.exit(0)

with open(infile) as file:
    target = os.path.abspath(file.read())

if not os.path.isfile(target):
    print 'Target file not found:\n{0}'.format(os.path.abspath(target))
    sys.exit(0)

print target
parent_dir = os.path.dirname(target)
name = getname(target)

record = TCSeqRec(target, file_format='phylip', name=name)
record.get_phyml_tree(model=model,ncat=ncat,datatype=datatype)

cPickle.dump(record, open('{0}/{1}.pickle'.format(parent_dir, name),'w'))
