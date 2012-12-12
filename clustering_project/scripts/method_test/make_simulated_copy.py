#!/usr/bin/env python

################################################################################
# Read in a SequenceCollection from disk and dump a simulated copy
################################################################################

################################################################################
# Commandline args:
#       -i      = input filename
#       -o      = output filename
#       -t      = temp dir
################################################################################

import argparse
import re
import sys
from errors import filecheck_and_quit, directorycheck_and_make

progname    = re.compile('[A-Za-z0-9.-_]+').search(sys.argv[0]).group()
desc        = 'Read in a SequenceCollection from disk and dump a simulated copy'
input_help  = 'Filepath+name of gzipped SequenceCollection object'
output_help = 'Output path+filename'

parser = argparse.ArgumentParser(prog=progname, description=desc)
parser.add_argument('-i', dest='input_file', help=input_help, type=str)
parser.add_argument('-o', dest='output_file', help=output_help, type=str)
parser.add_argument('-t', dest='tmpdir', type=str, default=None)

args        = parser.parse_args()
input_file  = args.input_file
output_file = args.output_file
tmpdir      = args.tmpdir

filecheck_and_quit(input_file)
directorycheck_and_make(output_dir)

from sequence_collection import SequenceCollection

sc = SequenceCollection.gunzip(input_file)
