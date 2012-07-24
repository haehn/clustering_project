#!/usr/bin/python
# -*- coding: utf-8 -*-

import re
import os
import dendropy as dpy
import numpy as np
import ete2
import random
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
        self.output = output

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

    def pam2sps(self, multiplier=0.01):
        """
        Scales branch lengths by an order of `multiplier`.
        Default is 0.01, converting PAM units to substitutions 
        per site.
        multiplier = 'sps2pam' scales by 100, performing the 
        opposite operation.
        multiplier = 'strip' removes branch lengths entirely
        """

        reg_ex = re.compile('(?<=:)[0-9.]+')

        converter = lambda a: str(multiplier * float(a.group()))
        strip_lengths = lambda d: ''

        input_string = self.newick

        if multiplier == 'pam2sps':
            multiplier = 0.01
        elif multiplier == 'sps2pam':
            multiplier = 100

        # Set the output string according to selection

        if multiplier == 'strip':
            output_string = reg_ex.sub(strip_lengths,
                    input_string).replace(':', '')
        else:

            output_string = reg_ex.sub(converter, input_string)

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
            if not self.newick.endswith('\n'):
                self.newick += '\n'
            writer.write(self.newick)
        writer.close()
        return outfile

    def run_phyml(
        self,
        model,
        alignment_file,
        datatype,
        name,
        ncat=4,
        verbose=True,
        overwrite=True,
        ):

        if not overwrite and self.newick:
            return self
        command = \
            'phyml -m {0} -i {1} -d {2} -c {3} -a e -b 0 --sequential > /dev/null'.format(model,
                alignment_file, datatype, ncat)
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

    def run_bionj(self, model, alignment_file, datatype, ncat, name, overwrite=True):

        if not overwrite and self.newick:
            return self
        command = 'phyml -m {0} -i {1} -d {2} -c {3} -b 0 -o n --sequential > /dev/null'.format(model,
                      alignment_file, datatype, ncat)
        process = Popen(command, shell=True, stdin=PIPE, stdout=PIPE,
                        stderr=PIPE).wait()
        tree = open('{0}_phyml_tree.txt'.format(alignment_file)).read()
        score = float(re.compile('(?<=Log-likelihood: ).+'
                      ).search(open('{0}_phyml_stats.txt'.format(alignment_file)).read()).group())
        output = \
            open('{0}_phyml_stats.txt'.format(alignment_file)).read()
        os.system('rm {0}_phyml_tree.txt {0}_phyml_stats.txt'.format(alignment_file))  # Cleanup

        (self.newick, self.score, self.program, self.name,
         self.output) = (tree, score, 'bionj', name, output)
        return Tree(tree, score, 'bionj', name, output)

    def run_raxml(
        self,
        model,
        alignment_file,
        name,
        tmpdir,
        guide=False,
        overwrite=True,
        ):

        if not overwrite and self.newick:
            return self
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
            output = None
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
        overwrite=True,
        ):

        if not overwrite and self.newick:
            return self
        command = 'TreeCollection {0} {1} {2} {3}'.format(dv_file,
                map_file, label_file, tree_file)
        process = Popen(command, shell=True, stdin=PIPE, stdout=PIPE,
                        stderr=PIPE)
        (stdout, stderr) = process.communicate()
        info = stdout.split()
        tree = info[-2]
        score = float(info[-1])

        (self.newick, self.score, self.program, self.name,
         self.output) = (tree, score, 'TreeCollection', name, stdout)
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

    def random_topology(
        self,
        nspecies,
        names=None,
        rooted=False,
        ):
        """
        Use ete2 to make a random topology
        Then add random branch lengths drawn from
        some distribution (default = gamma)
        Inner and leaf edge lengths can be drawn from differently parameterised
        versions of the distribution
        """

        if names:
            random.shuffle(names)
        t = ete2.Tree()
        t.populate(nspecies, names_library=names)
        if rooted:
            t.set_outgroup(t.children[0])
        else:
            t.unroot()

        t_as_newick = t.write()
        t_as_newick = t_as_newick.replace(')1', ')')
        return Tree(t_as_newick, name='random tree').pam2sps('strip')

    def randomise_branch_lengths(
        self,
        inner_edges,
        leaves,
        distribution_func=np.random.gamma,
        output_format=5
        ):
        """
        inner_edges and leaves are tuples describing the parameters of the
        distribution function
        distribution_func is a function generating samples from a probability
        distribution (eg gamma, normal ...)
        """
        t = ete2.Tree(self.newick)

        for node in t.iter_descendants():
            if node.children:
                node.dist = max(0,distribution_func(*inner_edges)) # 0 length or greater
            else:
                node.dist = max(0,distribution_func(*leaves))

        t_as_newick = t.write(format=output_format)
        return Tree(t_as_newick, name='random tree')

    def scale(self, scaling_factor):

        return self.pam2sps(scaling_factor)

    def nni(self):

        t = ete2.Tree(self.newick)
        original_outgroup = t.children[0]

        # Pick edge randomly by choosing a node and its parent

        inner_nodes = list(set(t.iter_descendants())
                           - set(t.iter_leaves()))
        node1 = random.choice(inner_nodes)
        t.set_outgroup(node1)
        node2 = node1.get_sisters()[0]
        child1 = random.choice(node1.children)
        child2 = random.choice(node2.children)

        child1.detach()
        child2.detach()
        node1.add_child(child2)
        node2.add_child(child1)

        t.set_outgroup(original_outgroup)
        t_as_newick = t.write(format=5)

        return Tree(t_as_newick, self.score, self.program, self.name,
                    self.output)
