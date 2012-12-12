#!/usr/bin/env python

################################################################################
# Reads a SequenceCollection object from a pickle, makes a randomised copy, 
# dumps records to disk ready for phylogenetic inference
################################################################################

################################################################################
# Commandline args:
#      -i = input file: a python pickle serialization of
#                       a SequenceCollection object
#                       containing a set of sequences
#
#      -o = output directory: the randomised copies of the sequences
#                       will be dumped in here, in phylip format.
#                       Names are hashed, so a translation dictionary
#                       is pickled here, along with the randomised
#                       copy of the original SequenceCollection object
#
# Outputs:
#      <files>.phy
#      hash_translation.pkl
#      scrand.pkl
################################################################################


from errors import filecheck_and_quit, directorycheck_and_make
import cPickle
import sys
import argparse
import re

progname = re.compile('[A-Za-z0-9.-_]+').search(sys.argv[0]).group()

desc = '\n'.join(['Read a SequenceCollection from pickle,',
                 'make a randomised copy,', 'dump records'])
input_help = 'Path+Filename for the input pickle file'
output_help = \
    'Path to output directory. Will be created if doesn\'t exist'
parser = argparse.ArgumentParser(prog=progname, description=desc,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-i', '--input', help=input_help, type=str)
parser.add_argument('-o', '--output', help=output_help, type=str)

args = vars(parser.parse_args())

pickle = args['input']
output_dir = args['output']

filecheck_and_quit(pickle)           # can't find file -> quit
directorycheck_and_make(output_dir)  # can't find directory -> create it

sc = cPickle.load(file(pickle))

scrand = sc.make_randomised_copy()

scrand.dump_records(output_dir)

cPickle.dump(scrand, open('{0}/scrand.pkl'.format(output_dir), 'w'))
