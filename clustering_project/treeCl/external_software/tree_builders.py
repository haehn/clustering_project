#!/usr/bin/env python

from . import ExternalSoftware
from ..errors import FileError, filecheck_and_raise
from ..utils.file_utils import *

local_dir = path_to(__file__)


class Phyml(ExternalSoftware):

    default_binary = 'phyml'

    def call(self, verbose=False):
        cmd = ' '.join([self.binary] + ['{0} {1}'.format(k, v) for (k,
                       v) in self.flags.items()])
        (stdout, stderr) = subprocess(cmd)
        if verbose:
            print stdout, stderr
        return (stdout, stderr)

    def read(self, filename):
        tree_filename = filename + '_phyml_tree.txt'
        stats_filename = filename + '_phyml_stats.txt'
        self.add_tempfile(tree_filename)
        self.add_tempfile(stats_filename)
        with open(tree_filename) as treefile:
            with open(stats_filename) as statsfile:
                return (treefile.read(), statsfile.read())

    def run(self, sequence_record):
        filename = self.write(sequence_record)
        filecheck_and_raise(filename)
        self.add_flag('-i', filename)
        print 'Running phyml on', sequence_record.name
        (stdout, stderr) = self.call()
        (tree, stats) = self.read(filename)
        self.clean()
        return (tree, stats)

    def write(self, sequence_record):
        filename = ((sequence_record.name if len(sequence_record.name)
                    < 40 else sequence_record.hashname()) if sequence_record.name else 'tmp_phyml_input'
                    )

        filename = '{0}/{1}.phy'.format(self.tmpdir, filename)
        sequence_record.write_phylip(filename)
        self.add_tempfile(filename)
        return filename


class TreeCollection(ExternalSoftware):

    default_binary = 'TreeCollection'
