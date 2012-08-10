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
infile = args['infile'] + os.environ['LSB_JOBINDEX']
# infile = sys.argv[1]
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
if 'dna_alignments' in target:
    record.get_phyml_tree(model='GTR',ncat=4,datatype='nt')
elif 'aa_alignments' in target:
    record.get_phyml_tree(model='WAG',ncat=4,datatype='aa')
cPickle.dump(record, open('{0}/{1}.pickle'.format(parent_dir, name),'w'))
