#!/usr/bin/python
# -*- coding: utf-8 -*-

import re
import os
import dendropy as dpy
from subprocess import Popen, PIPE, call


class Tree(object):

    """
    Class for storing the results of phylogenetic inference
    """

    def __init__(
        self,
        newick=None,
        score=0,
        program=None,
        name=None,
        output=None,
        ):

        self.newick = newick
        self.score = score
        self.program = program
        self.name = name

    def __str__(self):
        """
        Represents the object's information inside
        a newick comment, so is still interpretable
        by a (good) newick parser
        """

        s = '[Tree Object:\n'
        if self.name:
            s += 'Name:\t' + self.name + '\n'
        s += 'Program:\t{0}\n'.format(self.program) \
            + 'Score:\t{0}\n'.format(self.score) \
            + 'Tree:\t]{0}\n'.format(self.newick)
        return s

    def pam2sps(self, conversion='sps2pam'):
        """
        Scales branch lengths by an order of 100, converting
        PAM units to substitutions per site, and vice versa
        """

        reg_ex = re.compile('(?<=:)[0-9.]+')

        convert_pam_to_sps = lambda a: str(0.01 * float(a.group()))
        convert_sps_to_pam = lambda b: str(100 * float(b.group()))
        strip_lengths = lambda c: ''

        input_string = self.newick
        if conversion == 'pam2sps':
            output_string = reg_ex.sub(convert_pam_to_sps, input_string)
        elif conversion == 'sps2pam':
            output_string = reg_ex.sub(convert_sps_to_pam, input_string)
        else:
            output_string = reg_ex.sub(strip_lengths,
                    input_string).replace(':', '')

        return Tree(output_string, self.score, self.program, self.name)

    def read_from_file(self, infile, name=None):
        """
        This and the write_to_file function allow the class to be
        easily stored and reconstituted without using a pickle or
        JSON
        """

        program = None
        tree = None
        score = None
        self.name = name
        reader = open(infile)
        try:
            for line in reader:
                line = [l.rstrip().replace(']', '') for l in
                        line.split()]
                if not name and line[0] == 'Name:':
                    self.name = line[1]
                elif line[0] == 'Program:':
                    self.program = line[1]
                elif line[0] == 'Tree:':
                    self.newick = line[1]
                elif line[0] == 'Score:':
                    self.score = line[1]
        except IndexError:
            return
        return self

    def write_to_file(self, outfile, metadata=False):
        """
        Writes a string representation of the object's contents
        to file. This can be read using read_from_file to 
        reconstruct the Tree object, if metadata is included (i.e.
        metadata=True)
        """

        writer = open(outfile, 'w')
        if metadata:
            writer.write(str(self))
        else:
            writer.write(self.newick)
        writer.close()

    def run_phyml(
        self,
        model,
        alignment_file,
        datatype,
        name,
        verbose=True,
        ):

        command = \
            'phyml -m {0} -i {1} -d {2} -c 4 -a e -b 0 --sequential > /dev/null'.format(model,
                alignment_file, datatype)
        if verbose:
            command = command.replace('> /dev/null', '')
        process = Popen(command, shell=True, stdin=PIPE, stdout=PIPE,
                        stderr=PIPE).wait()
        tree = open('{0}_phyml_tree.txt'.format(alignment_file)).read()
        score = float(re.compile('(?<=Log-likelihood: ).+'
                      ).search(open('{0}_phyml_stats.txt'.format(alignment_file)).read()).group())
        output = \
            open('{0}_phyml_stats.txt'.format(alignment_file)).read()
        os.system('rm {0}_phyml_tree.txt {0}_phyml_stats.txt'.format(alignment_file))  # Cleanup
        (self.newick, self.score, self.program, self.name,
         self.output) = (tree, score, 'phyml', name, output)
        return Tree(tree, score, 'phyml', name, output)

    def run_raxml(
        self,
        model,
        alignment_file,
        name,
        tmpdir,
        guide=False,
        ):

        command = 'raxml -m {0} -s {1} -n {2} -w {3}'.format(model,
                alignment_file, name, tmpdir)
        if guide:
            command += ' -y'
        process = call(command, shell=True, stdin=PIPE, stdout=PIPE,
                       stderr=PIPE)
        if guide:
            dpy_tree = dpy.Tree()
            dpy_tree.read_from_stream(open('{0}/RAxML_parsimonyTree.{1}'.format(tmpdir,
                    name)), 'newick')
            dpy_tree.resolve_polytomies()
            for n in dpy_tree.postorder_node_iter():
                n.edge_length = 1
            tree = dpy_tree.as_newick_string()
            if not tree.rstrip().endswith(';'):
                tree += ';\n'
            score = None
        else:
            tree = open('{0}/RAxML_bestTree.{1}'.format(tmpdir,
                        name)).read()
            score = float(re.compile('(?<=Score of best tree ).+'
                          ).search(open('{0}/RAxML_info.{1}'.format(tmpdir,
                          name)).read()).group())
            output = open('{0}/RAxML_info.{1}'.format(tmpdir,
                          name)).read()
        os.system('rm {0}/*.{1}'.format(tmpdir, name))  # Cleanup
        (self.newick, self.score, self.program, self.name,
         self.output) = (tree, score, 'raxml', name, output)
        return Tree(tree, score, 'raxml', name, output)

    def run_treecollection(
        self,
        dv_file,
        map_file,
        label_file,
        tree_file,
        name,
        ):

        command = 'TreeCollection {0} {1} {2} {3}'.format(dv_file,
                map_file, label_file, tree_file)
        process = Popen(command, shell=True, stdin=PIPE, stdout=PIPE,
                        stderr=PIPE)
        (stdout, stderr) = process.communicate()
        info = stdout.split()
        tree = info[-2]
        score = float(info[-1])
        (self.newick, self.score, self.program, self.name) = (tree,
                score, 'TreeCollection', name, stdout)
        return Tree(tree, score, 'TreeCollection', name,
                    stdout).pam2sps('pam2sps')

    def _unpack_raxml_args(packed_args):
        """
        Used as helper to enable parallelising multiple system calls
        to raxml
        """

        return run_raxml(*packed_args)  # * to unpack

    def _unpack_phyml_args(packed_args):
        """
        Used as helper to enable parallelising multiple system calls
        to phyml
        """

        return run_phyml(*packed_args)  # * to unpack

    def _unpack_TC_args(packed_args):
        """
        Used as helper to enable parallelising multiple system calls
        to treecollection
        """

        return run_treecollection(*packed_args)  # * to unpack
