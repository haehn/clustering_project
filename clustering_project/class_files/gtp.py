#!/usr/bin/env python

# from sequence_record import TCSeqRec


from errors import FileError
import numpy as np
import os

local_dir = os.path.dirname(__file__)

class GTP(object):

    """
    Interacts with gtp.jar to calculate geodesic distances from trees
    """

    def __init__(self, binary='', tmpdir='/tmp'):

        supplied_binary = binary.rstrip('/')
        
        if self.can_locate(supplied_binary):
            self.binary = supplied_binary
        
        else:
            default_binary = self.locate_binary()
            self.binary = default_binary

        if self.binary is None:
            raise FileError(supplied_binary)

        self.tmpdir = tmpdir.rstrip('/')

    def can_locate(self, binary):
        return os.path.isfile(binary)

    def locate_binary(self):
        if 'GTP_PATH' in os.environ:
            print 'Using $GTP_PATH to find gtp.jar:',
            binary = '/'.join((os.environ['GTP_PATH'].rstrip('/'),
                              'gtp.jar'))
        else:
            print 'Using this as gtp.jar:',
            binary = '/'.join((local_dir, 'gtp.jar'))
        print binary

        return (binary if self.can_locate(binary) else None)

    def call(self, rooted):

        bincall = 'java -jar {0}'.format(self.binary)
        if not rooted:
            flags = '-u -o'
        else:
            flags = '-o'
        outf = '{0}/output.txt'.format(self.tmpdir)

        inf = '{0}/geotrees.nwk'.format(self.tmpdir)

        cmd = ' '.join((bincall, flags, outf, inf))

        os.system(cmd)

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

    def clean(self):
        os.remove('{0}/geotrees.nwk'.format(self.tmpdir))
        os.remove('{0}/output.txt'.format(self.tmpdir))

    def writetmp(self, trees):
        with open('{0}/geotrees.nwk'.format(self.tmpdir), 'w') as tmpf:
            tmpf.write('\n'.join(tree.newick.rstrip() for tree in
                       trees))

    def allrooted(self, trees):
        return all(tree.rooted for tree in trees)

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

    def pairwise(self, tree1, tree2):
        return self.run( (tree1, tree2))[0,1]
        