#!/usr/bin/env python

################################################################################
# Read in a SequenceCollection from disk and dump records
################################################################################

################################################################################
# Commandline args:
#       -i      = input filename
#       -o      = output directory
#       -c      = dump post-clustering concatenated records
#               instead of pre-clustering single records
################################################################################

import argparse
import re
import sys
from errors import filecheck_and_quit, directorycheck_and_make

progname    = re.compile('[A-Za-z0-9.-_]+').search(sys.argv[0]).group()
desc        = 'Read in a SequenceCollection from disk and dump records'
input_help  = 'Filepath+name of gzipped SequenceCollection object'
output_help = 'Directory to dump files in'
choice_help = \
    '\n'.join(['Choose to dump post-clustering concatenated records',
              'instead of pre-clustering single records'])
parser = argparse.ArgumentParser(prog=progname, description=desc)
parser.add_argument('-i', dest='input_file', help=input_help, type=str)
parser.add_argument('-o', dest='output_dir', help=output_help, type=str)
parser.add_argument('-c', dest='cluster_recs', action='store_true')

args         = parser.parse_args()
input_file   = args.input_file
output_dir   = args.output_dir.rstrip('/')
cluster_recs = args.cluster_recs

filecheck_and_quit(input_file)
directorycheck_and_make(output_dir)

from sequence_collection import SequenceCollection

sc = SequenceCollection.gunzip(input_file)
if cluster_recs:
    records = sc.get_cluster_records()
    sc.dump_records(output_dir, records)
else:
    sc.dump_records(output_dir)
