#!/usr/bin/env python

import argparse
import os
import sys


def fpath(s):
    """
    Helper function used when passing filepath arguments with argparse module.
    Trims all '/' characters from the end of the path string.
    """

    while s.endswith('/'):
        s = s[:-1]
    return s


desc = 'Splits TreeCollection DistVar.txt file into individual files'
parser = argparse.ArgumentParser(prog='distvar_splitter.py',
                                 description=desc,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-d', '--distvar',
                    help='Location of input distance-variance file',
                    type=fpath)
parser.add_argument('-g', '--genemap', 
                    help='Location of input genemap file',
                    type=fpath)
parser.add_argument('-p', '--prefix',
                    help='Optional prefix for output files', type=str,
                    default='gene_')

# parser.add_argument('-o', '--outdir')

args = vars(parser.parse_args())
dvfile = args['distvar']
gmfile = args['genemap']
prefix = args['prefix']

# Check infiles are real

for infile in (dvfile, gmfile):
    if not os.path.isfile(infile):
        print 'Can\'t find this file: {0}'.format(infile)
        sys.exit(0)


def get_dv_reader(f):
    reader = open(f)
    num_records = int(reader.readline())
    return (reader, num_records)

def get_gm_reader(f):
    reader = open(f)
    num_records, num_taxa = (int(x) for x in reader.readline().rstrip().split())
    return (reader, num_records, num_taxa)

def read_a_gm_record(reader, num_records, num_taxa):
    header = '1 {0}\n'.format(num_taxa)
    gm = reader.readline()
    if len(gm.split()) < num_taxa:
        print 'Line is of unexpected length'
        print 'Line: {0}'.format(gm)
        print 'Expected {0} taxa'.format(num_taxa)
        print 'Actual number of taxa: {0}'.format(len(gm.split()))
        sys.exit()
    return (header + gm)   

def read_a_dv_record(reader, num_records):
    (i, j, k) = (int(x) for x in reader.readline().split())
    header = '''1
{0} {1} 1
'''.format(i, j)
    matrix = []
    for _ in range(i):
        line = reader.readline().rstrip()
        if len(line.split()) != j:
            print 'Line is of unexpected length'
            print 'Line: {0}'.format(line)
            print 'Expected length = {0}'.format(len(line.split()))
            print 'Actual length = {0}'.format(j)
            sys.exit()
        matrix.append(line)
    matrix = '\n'.join(matrix) + '\n'
    if k == num_records:
        reader.close()
    return (header + matrix, k)


if __name__ == '__main__':
    (r, n) = get_dv_reader(dvfile)
    (gr, gn, t) = get_gm_reader(gmfile)
    if not n == gn:
        print 'distance-variance file and genemap file'
        print 'contain different numbers of sequences'
        sys.exit()
    while not r.closed:
        (matrix, number) = read_a_dv_record(r, n)
        name = '{0}{1:0>3}'.format(prefix,number)
        gm = read_a_gm_record(gr, gn, t)
        with open('{0}.dv'.format(name),'w') as dvwriter:
            dvwriter.write(matrix)
        with open('{0}.gm'.format(name),'w') as gmwriter:
            gmwriter.write(gm)

