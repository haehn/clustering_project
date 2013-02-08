#!/usr/bin/env python

from errors import FileError, filecheck_and_raise
from file_utils import *
import numpy as np
import hashlib

local_dir = path_to(__file__)


def hashname(name):
    return hashlib.sha1(name).hexdigest()


class ExternalSoftware(object):

    flags = {}
    default_binary = ''
    default_env = ''

    def __init__(self, supplied_binary='', tmpdir='/tmp'):

        if can_locate(supplied_binary):
            self.binary = supplied_binary
        else:

            default_binary = locate_file(self.default_binary,
                    self.default_env, local_dir)
            self.binary = default_binary

        if self.binary is None:
            raise FileError(supplied_binary)

        self.tmpdir = tmpdir.rstrip('/')

    def __str__(self):
        desc = 'Wrapper for {0}'.format(self.default_binary)
        return desc

    def add_flag(self, flag, value):
        self.flags[flag] = value

    def remove_flag(self, flag):
        del self.flags[flag]

    def call(self):
        pass

    def clean(self):
        pass

    def read(self):
        pass

    def run(self):
        pass

    def writetmp(self):
        pass


class GTP(ExternalSoftware):

    """
    Interacts with gtp.jar to calculate geodesic distances from trees
    """

    default_binary = 'gtp.jar'
    default_env = 'GTP_PATH'

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


class Phyml(ExternalSoftware):

    default_binary = 'phyml'

    def clean(self, filename):
        prefix = '{0}/{1}.phy_phyml_'.format(self.tmpdir, filename)
        delete(prefix + 'tree.txt')
        delete(prefix + 'stats.txt')

    def call(self):
        cmd = ' '.join([self.binary] + ['{0} {1}'.format(k,
                                v) for k,v in self.flags.items()])
        return subprocess(cmd)

    def writetmp(self, sequence_record):
        if '-i' in self.flags:
            filename = (self.flags['-i'] if len(self.flags['-i'])
                        < 100 else hashname(self.flags['-i']))
        else:
            print 'No input file to write!'
            return

        sequence_record.write_phylip('{0}/{1}.phy'.format(self.tmpdir,
                filename))
        return filename

    def read(self, filename):
        prefix = '{0}/{1}.phy_phyml_'.format(self.tmpdir, filename)
        with open(prefix + 'tree.txt') as treefile:
            with open(prefix + 'stats.txt') as statsfile:
                return (treefile.read(), statsfile.read())

    def run(self, sequence_record):
        filename = self.writetmp(sequence_record)
        filecheck_and_raise('{0}/{1}.phy'.format(self.tmpdir,
                filename))
        print 'Running phyml on', sequence_record.name
        self.call()


class TreeCollection(ExternalSoftware):

    default_binary = 'TreeCollection'
