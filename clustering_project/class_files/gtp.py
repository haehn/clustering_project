#!/usr/bin/env python

# from sequence_record import TCSeqRec

from errors import FileError
from file_utils import *
import numpy as np

local_dir = path_to(__file__)


class GTP(object):

    """
    Interacts with gtp.jar to calculate geodesic distances from trees
    """

    def __init__(self, supplied_binary='', tmpdir='/tmp'):

        if can_locate(supplied_binary):
            self.binary = supplied_binary
        else:

            default_binary = locate_file('gtp.jar', 'GTP_PATH',
                    local_dir)
            self.binary = default_binary

        if self.binary is None:
            raise FileError(supplied_binary)

        self.tmpdir = tmpdir.rstrip('/')

    def __str__(self):
        desc = 'Wrapper for gtp.jar - geodesic distance calculator'
        authors = '(Owen, Megan, and J Scott Provan. 2011.'
        title = \
            'A Fast Algorithm for Computing Geodesic Distances in Tree Space,'
        doi = 'doi:10.1109/TCBB.2010.3)'

        details = \
            'Jarfile: {0}\nTemp directory: {1}'.format(self.binary,
                self.tmpdir)

        return '\n'.join((
            desc,
            authors,
            title,
            doi,
            '',
            details,
            '',
            ))

    def allrooted(self, trees):
        return all(tree.rooted for tree in trees)

    def call(self, rooted):

        bincall = 'java -jar {0}'.format(self.binary)
        if not rooted:
            flags = '-u -o'
        else:
            flags = '-o'
        outf = '{0}/output.txt'.format(self.tmpdir)
        inf = '{0}/geotrees.nwk > /dev/null'.format(self.tmpdir)

        cmd = ' '.join((bincall, flags, outf, inf))
        syscall(cmd)

    def clean(self):
        delete('{0}/geotrees.nwk'.format(self.tmpdir))
        delete('{0}/output.txt'.format(self.tmpdir))

    def pairwise(self, tree1, tree2):
        return self.run((tree1, tree2))[0, 1]

    def read(self, size):
        matrix = np.zeros((size, size))

        try:
            with open('{0}/output.txt'.format(self.tmpdir)) as outf:
                for line in outf:
                    line = line.rstrip()
                    if line:
                        (i, j, value) = line.split()
                        i = int(i)
                        j = int(j)
                        value = float(value)
                        matrix[i, j] = matrix[j, i] = value

            return matrix
        except IOError, e:

            print 'There was an IOError: {0}'.format(e)
            print 'Geodesic distances couldn\'t be calculated'
            raise

    def run(self, trees):
        self.writetmp(trees)
        rooted = self.allrooted(trees)
        self.call(rooted)
        try:
            matrix = self.read(len(trees))
            self.clean()
            return matrix
        except IOError:
            print 'except'
            matrix = None
            raise

    def writetmp(self, trees):
        with open('{0}/geotrees.nwk'.format(self.tmpdir), 'w') as tmpf:
            tmpf.write('\n'.join(tree.newick.rstrip() for tree in
                       trees))


if __name__ == '__main__':
    from tree import Tree
    trees = [Tree.new_random_coal(10) for _ in range(100)]
    g = GTP()
    print g
    m = g.run(trees)
    print m
