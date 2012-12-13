#!/usr/bin/env python

################################################################################
# Read in a SequenceCollection from disk and print scores
################################################################################

################################################################################
# Commandline args:
#       -i      = input filename
#       -t      = phyml tree results directory
################################################################################

import argparse
import re
import sys
from errors import filecheck_and_quit, directorycheck_and_quit

progname = re.compile('[A-Za-z0-9.-_]+').search(sys.argv[0]).group()
desc = 'Read in a SequenceCollection from disk and print scores'
input_help = 'Filepath+name of gzipped SequenceCollection object'

parser = argparse.ArgumentParser(prog=progname, description=desc)
parser.add_argument('-i', dest='input_file', help=input_help, type=str)
parser.add_argument('-t', dest='phyml_dir', help=input_help, type=str)

args = parser.parse_args()
input_file = args.input_file
phyml_dir = args.phyml_dir.rstrip('/')

filecheck_and_quit(input_file)
directorycheck_and_quit(phyml_dir)

from sequence_collection import SequenceCollection

sc = SequenceCollection.gunzip(input_file)
cluster_records = sc.get_cluster_records()
sc.load_phyml_results(phyml_dir, records=cluster_records, use_hashname=True)
sc.update_scores()
sc.gzip(input_file)
scores = sorted(sc.get_scores(), key=lambda x: (x[0]))
for score in scores:
    print score