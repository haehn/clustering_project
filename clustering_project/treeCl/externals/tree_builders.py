#!/usr/bin/env python

from . import ExternalSoftware
from ..errors import FileError, filecheck_and_raise
from ..utils import dpy, fileIO
from treeCl.tree import Tree, Manipulations
import re

local_dir = fileIO.path_to(__file__)

class Result(object):

    def __init__(self, output, score, tree, name, program):
        self.output  = output
        self.score   = score
        self.tree    = tree
        self.name    = name
        self.program = program

    def __str__(self):
        return self.output

class Phyml(ExternalSoftware):

    default_binary = 'phyml'   
    score_regex = re.compile('(?<=Log-likelihood: ).+') 

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

        score = float(self.score_regex.search(stats).group(0))
        self.clean()
        return Result(stats, score, tree, sequence_record.name, 
            fileIO.basename(self.binary))

    def write(self, sequence_record):
        filename = sequence_record.get_name(default='tmp_phyml_input')
        filename = '{0}/{1}.phy'.format(self.tmpdir, filename)
        sequence_record.write_phylip(filename)
        self.add_tempfile(filename)
        return filename

    def get_datatype(self, sequence_record, datatype=None):
        datatype = datatype or sequence_record.datatype
        if datatype == 'protein':
            return {'-d': 'aa', '-m': 'WAG'}
        elif datatype == 'dna':
            return {'-d': 'nt', '-m': 'GTR'}

    def set_default_flags(
        self,
        sequence_record,
        analysis='ml',
        datatype=None,
        ):

        defaults = self.get_datatype(sequence_record, datatype=datatype)
        if defaults:
            defaults['-a'] = 'e'
            defaults['-b'] = 0
            defaults['-c'] = 4
            defaults['-q'] = ''
            if analysis == 'ml' or analysis == 'full':
                defaults['-o'] = 'tlr'
            elif analysis == 'nj' or analysis == 'bionj':
                defaults['-o'] = 'n'

            for flag in defaults:
                self.add_flag(flag, defaults[flag])


class TreeCollection(ExternalSoftware):

    default_binary = 'TreeCollection'

    def parse(self, output):
        output = output.split()
        score = float(output[-1])
        tree = output[-2]
        return (score, tree)

    def run(
        self,
        sequence_record,
        guidetree=None,
        datatype='protein',
        verbose=False,
        ):
        guidetree = guidetree or self.write_guidetree(sequence_record, datatype)
        tmpfiles = self.write(sequence_record)
        self.add_flag('-D', tmpfiles['dv'])
        self.add_flag('-M', tmpfiles['map'])
        self.add_flag('-L', tmpfiles['lab'])
        self.add_flag('-T', guidetree)
        (stdout, stderr) = self.call()
        self.clean()
        if verbose: 
            print stdout, stderr
        score, tree = self.parse(stdout)
        tree = Manipulations(tree).pam2sps()
        return Tree(tree, score, fileIO.basename(self.binary),
                sequence_record.name, stdout)

    def write(self, sequence_record):
        """ Write the distance-variance (dv) file, the labels file and the map
        file all together, because they share information, and this way we only
        have to look it up once"""

        # Look up info

        dv_info = sequence_record.dv
        num_matrices = len(dv_info)
        if num_matrices == 0:
            print 'No distance-variance matrix available'
            return
        all_labels = sequence_record.headers
        labels_len = len(all_labels)

        base_filename = sequence_record.get_name(default='tmp_TC')
        f = {}
        f['dv']  = '{0}/{1}_dv.txt'.format(self.tmpdir, base_filename)
        f['lab'] = '{0}/{1}_labels.txt'.format(self.tmpdir, base_filename)
        f['map'] = '{0}/{1}_map.txt'.format(self.tmpdir, base_filename)
        dv_file     = open(f['dv'], 'w')
        labels_file = open(f['lab'], 'w')
        map_file    = open(f['map'], 'w')

        # Write headers
        dv_file.write('{0}\n'.format(num_matrices))
        map_file.write('{0} {1}\n'.format(num_matrices, labels_len))
        labels_file.write('{0}\n{1}\n'.format(labels_len, ' '.join(all_labels)))
        labels_file.close()

        # Write info
        for (i, (matrix, labels)) in enumerate(dv_info, start=1):
            labels = labels.split()
            dim = len(labels)
            dv_file.write('{0} {0} {1}\n{2}\n'.format(dim, i, matrix))
            dv_file.flush()
            for lab in all_labels:
                map_file.write(('{0} '.format(labels.index(lab) + 1) if lab
                               in labels else '-1 '))
            map_file.write('\n')
            map_file.flush()
        dv_file.close()
        map_file.close()

        for k in f:
            self.add_tempfile(f[k]) # for cleanup

        return f

    def write_guidetree(self, sequence_record, datatype='protein'):
        base_filename = sequence_record.get_name(default='tmp_TC')
        bionj = Phyml()
        bionj.set_default_flags(sequence_record, analysis='bionj',
                                datatype=datatype)
        tree = bionj.run(sequence_record).tree
        tree = dpy.bifurcate_base(tree)
        filename = '{0}/{1}_tree.nwk'.format(self.tmpdir, base_filename)
        with open(filename, 'w') as tree_file:
            tree_file.write(tree)
        self.add_tempfile(filename)
        return filename
