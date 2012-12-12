#!/usr/bin/env python

################################################################################
# Submits jobs to LSF
################################################################################

################################################################################
# Commandline args:
#      -i   = input directory: contains sequence alignments in PHYLIP format
#                       (can be sequential or interleaved)
#
#      -d   = datatype  either 'protein' for protein data or 'dna' for dna
#                       data. (Sets phyml model to either WAG or GTR)
#
#      -m   = method:   either 'bionj' or 'phyml'
#
#      -int = interleaved: add this flag if the data are interleaved
#
# Outputs:
#      {<file>_phyml_tree.txt }  for each sequence file in the input directory,
#      {<file>_phyml_stats.txt}  written to the input directory
#
################################################################################

from errors import filecheck_and_quit, directorycheck_and_raise
import sys
import argparse
import glob
import os
import re

progname = re.compile('[A-Za-z0-9.-_]+').search(sys.argv[0]).group()

desc = 'Send off LSF jobs to get phyml trees for a set of records'
input_help = 'Path to input directory'
datatype_help = 'Datatype = \'protein\' || \'dna\''
method_help = 'Choose either \'phyml\' || \'bionj\''
seq_help = 'Sequences are interleaved'
datatype_choices = ['protein', 'dna']
method_choices = ['phyml', 'bionj']
formatter = argparse.ArgumentDefaultsHelpFormatter
parser = argparse.ArgumentParser(prog=progname, description=desc,
                                 formatter_class=formatter)
parser.add_argument('-i', '--input', help=input_help, type=str)
parser.add_argument(
    '-d',
    '--datatype',
    help=datatype_help,
    type=str,
    default='protein',
    choices=datatype_choices,
    )
parser.add_argument(
    '-m',
    '--method',
    help=method_help,
    type=str,
    default='phyml',
    choices=method_choices,
    )
parser.add_argument('-int', '--interleaved', help=seq_help,
                    action='store_true', default=False)

args = vars(parser.parse_args())
input_dir = args['input'].rstrip('/')
datatype = args['datatype']
method = args['method']
interleaved = args['interleaved']

directorycheck_and_raise(input_dir)

files = glob.glob('{0}/*.phy'.format(input_dir))

for f in files:
    if method == 'phyml':
        if datatype == 'protein':
            phyml_command = ' '.join(['phyml -m WAG -i {0}'.format(f),
                    '-d aa -c 4 -a e -b 0', '--no_memory_check'])
        elif datatype == 'dna':
            phyml_command = ' '.join(['phyml -m GTR -i {0}'.format(f),
                    '-d nt -c 4 -a e -b 0', '--no_memory_check'])
    elif method == 'bionj':

        if datatype == 'protein':
            phyml_command = ' '.join(['phyml -m WAG -i {0}'.format(f),
                    '-d aa -c 4 -b 0 -o n --no_memory_check'])
        elif datatype == 'dna':
            phyml_command = ' '.join(['phyml -m GTR -i {0}'.format(f),
                    '-d nt -c 4 -b 0 -o n --no_memory_check'])

    if not interleaved:
        phyml_command += ' --sequential'

    os.system('bsub -o /dev/null ' + phyml_command)
