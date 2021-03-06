#!/usr/bin/env python

################################################################################
# Read in a SequenceCollection from disk and make a randomised copy,
# then dump the alignments into an output directory
################################################################################

################################################################################
# Commandline args:
#       -i      = input filename (gzipped pickle [.pkl.gz])
#       -o      = output directory
#       -t      = temp dir
################################################################################

import argparse
import re
import sys
from errors import filecheck_and_quit, directorycheck_and_make

progname = re.compile('[A-Za-z0-9.-_]+').search(sys.argv[0]).group()
desc = \
    'Read in a SequenceCollection from disk and dump a randomised copy'
input_help = 'Filepath+name of gzipped SequenceCollection object'
output_help = 'Output directory'

parser = argparse.ArgumentParser(prog=progname, description=desc)
parser.add_argument('-i', dest='input_file', help=input_help, type=str)
parser.add_argument('-o', dest='output_dir', help=output_help, type=str)
parser.add_argument('-t', dest='tmpdir', type=str, default=None)

args = parser.parse_args()
input_file = args.input_file
output_dir = args.output_dir
tmpdir = args.tmpdir

filecheck_and_quit(input_file)
directorycheck_and_make(output_dir)

from sequence_collection import SequenceCollection

sc = SequenceCollection.gunzip(input_file)
scrand = sc.make_randomised_copy(tmpdir)
records = scrand.get_records()
scrand.dump_records(output_dir)
