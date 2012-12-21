#!/usr/bin/env python

################################################################################
# Convert raw scores to suitable dataframe, ready for ggplot
################################################################################

################################################################################
# Commandline args:
#       -i      = input filename
#       -o      = output filename
#       -c      = category (optional)
#       -n      = number (optional)
#       -header = header (true/false)
################################################################################

import argparse
import re
import sys
from errors import filecheck_and_quit, directorycheck_and_quit

progname = re.compile('[A-Za-z0-9.-_]+').search(sys.argv[0]).group()
desc = 'Convert raw scores to suitable dataframe, ready for ggplot'
category_choices = ['Observed', 'Randomised', 'Simulated', 'NA']

parser = argparse.ArgumentParser(prog=progname, description=desc)
parser.add_argument('-i', dest='input_file', type=str)
parser.add_argument('-o', dest='output_file', type=str)
parser.add_argument('-c', dest='category', choices=category_choices,
                    type=str, default='NA')
parser.add_argument('-n', dest='number', type=str, default='NA')
parser.add_argument('-header', dest='header', action='store_true',
                        default=False)

args = parser.parse_args()
input_file = args.input_file
output_file = args.output_file
category = args.category
number = str(args.number)
header = args.header

filecheck_and_quit(input_file)

from sequence_collection import SequenceCollection

sc = SequenceCollection.gunzip(input_file)
scores = sorted(sc.get_scores(), key=lambda x: x[0])

with open(output_file, 'w') as outf:
    if header:
        outf.write('\t'.join([
            'Cluster',
            'Dist',
            'Method',
            'LogLikelihood',
            'Diff',
            'Datatype',
            'Dataset',
            ]))
        outf.write('\n')

    prev_lk = None
    prev_meth = None
    for score in scores:
        ((dist, meth, clst), lk) = score
        if prev_meth != meth:
            prev_lk = None
            prev_meth = meth
        if prev_lk:
            df = lk - prev_lk
            prev_lk = lk
        else:
            df = 'NA'
            prev_lk = lk
        outf.write('\t'.join([str(clst), dist, meth,
                              str(lk), str(df), category,
                              number]))
        outf.write('\n')
