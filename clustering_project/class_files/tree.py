#!/usr/bin/env python

import re
import os
import dendropy as dpy
from dendropy import treesim
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
        rooted=None,
        ):

        self.newick = newick
        self.score = score
        self.program = program
        self.name = name
        self.output = output
        self.rooted = rooted

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
            + 'Rooted:\t{0}\n'.format(self.rooted) \
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

        return Tree(output_string, self.score, self.program, self.name,
                    self.rooted)

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

            writeable = self.newick
            if not writeable.endswith('\n'):
                writeable += '\n'
            writer.write(writeable)
        writer.close()
        return outfile

    @classmethod
    def check_rooted(cls, newick):
        t = dpy.Tree()
        t.read_from_string(newick, 'newick')
        root_degree = len(t.seed_node.child_nodes())
        return root_degree == 2

    @classmethod
    def deroot_tree(cls, newick):
        t = dpy.Tree()
        t.read_from_string(newick, 'newick')
        t.deroot()
        return t.as_newick_string() + ';'

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
            'phyml -m {0} -i {1} -d {2} -c {3} -a e -b 0 --sequential --no_memory_check #> /dev/null'.format(model,
                alignment_file, datatype, ncat)
        if verbose:
            command = command.replace('> /dev/null', '')
        process = Popen(command, shell=True, stdin=PIPE, stdout=PIPE,
                        stderr=PIPE)
        process.wait()
        tree = open('{0}_phyml_tree.txt'.format(alignment_file)).read()
        score = float(re.compile('(?<=Log-likelihood: ).+'
                      ).search(open('{0}_phyml_stats.txt'.format(alignment_file)).read()).group())
        output = \
            open('{0}_phyml_stats.txt'.format(alignment_file)).read()
        rooted = self.check_rooted(tree)
        os.system('rm {0}_phyml_tree.txt {0}_phyml_stats.txt'.format(alignment_file))  # Cleanup

        (
            self.newick,
            self.score,
            self.program,
            self.name,
            self.output,
            self.rooted,
            ) = (
            tree,
            score,
            'phyml',
            name,
            output,
            rooted,
            )
        return Tree(
            tree,
            score,
            'phyml',
            name,
            output,
            rooted,
            )

    def run_bionj(
        self,
        model,
        alignment_file,
        datatype,
        ncat,
        name,
        overwrite=True,
        ):

        if not overwrite and self.newick:
            return self
        command = \
            'phyml -m {0} -i {1} -d {2} -c {3} -b 0 -o n --sequential --no_memory_check #> /dev/null'.format(model,
                alignment_file, datatype, ncat)
        process = Popen(command, shell=True, stdin=PIPE, stdout=PIPE,
                        stderr=PIPE)
        process.wait()
        tree = open('{0}_phyml_tree.txt'.format(alignment_file)).read()
        score = float(re.compile('(?<=Log-likelihood: ).+'
                      ).search(open('{0}_phyml_stats.txt'.format(alignment_file)).read()).group())
        output = \
            open('{0}_phyml_stats.txt'.format(alignment_file)).read()
        rooted = self.check_rooted(tree)
        os.system('rm {0}_phyml_tree.txt {0}_phyml_stats.txt'.format(alignment_file))  # Cleanup

        (
            self.newick,
            self.score,
            self.program,
            self.name,
            self.output,
            self.rooted,
            ) = (
            tree,
            score,
            'bionj',
            name,
            output,
            rooted,
            )
        return Tree(
            tree,
            score,
            'bionj',
            name,
            output,
            rooted,
            )

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
        rooted = self.check_rooted(tree)
        os.system('rm {0}/*.{1}'.format(tmpdir, name))  # Cleanup

        (
            self.newick,
            self.score,
            self.program,
            self.name,
            self.output,
            self.rooted,
            ) = (
            tree,
            score,
            'raxml',
            name,
            output,
            rooted,
            )
        return Tree(
            tree,
            score,
            'raxml',
            name,
            output,
            rooted,
            )

    @classmethod
    def new_treecollection_tree(
        cls,
        dv_file,
        map_file,
        label_file,
        tree_file,
        name,
        overwrite=True,
        deroot=True,
        ):

        command = \
            'TreeCollection -D {0} -M {1} -L {2} -T {3}'.format(dv_file,
                map_file, label_file, tree_file)
        process = Popen(command, shell=True, stdin=PIPE, stdout=PIPE,
                        stderr=PIPE)
        process.wait()
        (stdout, stderr) = process.communicate()
        info = stdout.split()
        tree = info[-2]
        if deroot:
            tree = cls.deroot_tree(tree)
        score = float(info[-1])
        rooted = cls.check_rooted(tree)

        return Tree(
            tree,
            score,
            'TreeCollection',
            name,
            stdout,
            rooted,
            ).pam2sps('pam2sps')

    def run_treecollection(
        self,
        dv_file,
        map_file,
        label_file,
        tree_file,
        name,
        overwrite=True,
        deroot=True,
        ):

        if not overwrite and self.newick:
            return self

        new_tree = self.new_treecollection_tree(
            dv_file,
            map_file,
            label_file,
            tree_file,
            name,
            overwrite,
            deroot,
            )

        (
            self.newick,
            self.score,
            self.program,
            self.name,
            self.output,
            self.rooted,
            ) = (
            new_tree.newick,
            new_tree.score,
            new_tree.program,
            new_tree.name,
            new_tree.output,
            new_tree.rooted,
            )

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
        output_format=5,
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
                node.dist = max(0, distribution_func(*inner_edges))  # 0 length or greater
            else:
                node.dist = max(0, distribution_func(*leaves))

        t_as_newick = t.write(format=output_format)
        return Tree(t_as_newick, name='random tree')

    def scale(self, scaling_factor):

        return self.pam2sps(scaling_factor)

    def random_yule(self, nspecies=None, names=None):
        if names and nspecies:
            if not nspecies == len(names):
                nspecies = len(names)
        elif names and not nspecies:
            nspecies = len(names)
        elif not names and nspecies:
            names = ['Sp{0}'.format(i) for i in range(1, nspecies + 1)]
        elif not names and not nspecies:
            nspecies = 16
            names = ['Sp{0}'.format(i) for i in range(1, nspecies + 1)]

        taxon_set = dpy.TaxonSet(names)
        tree = treesim.uniform_pure_birth(taxon_set)

        newick = '[&R] ' + tree.as_newick_string()
        if not newick.endswith(';'):
            newick += ';'

        return Tree(newick)

    def random_coal(self, nspecies=None, names=None):
        if names and nspecies:
            if not nspecies == len(names):
                nspecies = len(names)
        elif names and not nspecies:
            nspecies = len(names)
        elif not names and nspecies:
            names = ['Sp{0}'.format(i) for i in range(1, nspecies + 1)]
        elif not names and not nspecies:
            nspecies = 16
            names = ['Sp{0}'.format(i) for i in range(1, nspecies + 1)]

        taxon_set = dpy.TaxonSet(names)
        tree = treesim.pure_kingman(taxon_set)

        newick = '[&R] ' + tree.as_newick_string()
        if not newick.endswith(';'):
            newick += ';'

        return Tree(newick)

    def get_constrained_gene_tree(
        self,
        scale_to=None,
        population_size=None,
        trim_names=True,
        ):
        """
        Using the current tree object as a species tree, generate
        a gene tree using the constrained Kingman coalescent
        process from dendropy.
        The species tree should probably be a valid, ultrametric
        tree, generated by some pure birth, birth-death or coalescent
        process, but no checks are made.
        Optional kwargs are: 
        -- scale_to, which is a floating point value to
        scale the total tree tip-to-root length to,
        -- population_size, which is a floating point value which 
        all branch lengths will be divided by to convert them to coalescent 
        units, and
        -- trim_names, boolean, defaults to true, trims off the number
        which dendropy appends to the sequence name
        """

        tree = dpy.Tree()
        tree.read_from_string(self.newick, 'newick')

        for leaf in tree.leaf_iter():
            leaf.num_genes = 1

        tree_height = tree.seed_node.distance_from_root() \
            + tree.seed_node.distance_from_tip()

        if scale_to:
            population_size = tree_height / scale_to

        for edge in tree.preorder_edge_iter():
            edge.pop_size = population_size

        gene_tree = treesim.constrained_kingman(tree)[0]

        if trim_names:
            for leaf in gene_tree.leaf_iter():
                leaf.taxon.label = leaf.taxon.label.replace('\'', ''
                        ).split('_')[0]

        newick = '[&R] ' + gene_tree.as_newick_string()
        if not newick.endswith(';'):
            newick += ';'

        return Tree(newick)

    def nni(self):

        # Function to perform a random nearest-neighbour interchange on a tree
        # using Dendropy

        # The dendropy representation of a tree is as if rooted (even when it's
        # not) The node which acts like the root is called the seed node, and this
        # can sit in the middle of an edge which would be a leaf edge in the
        # unrooted tree. NNI moves are only eligible on internal edges, so we
        # need to check if the seed node is sat on a real internal edge, or
        # a fake one.

        # Make a dendropy representation of the tree

        tree = dpy.Tree()
        tree.read_from_string(self.newick, 'newick')
        seed = tree.seed_node
        resolved = False
        if len(seed.child_nodes()) > 2:
            print 'Resolve root trisomy'
            tree.resolve_polytomies()
            resolved = True

        # Make a list of internal edges not including the root edge

        edge_list = list(tree.preorder_edge_iter(lambda edge: \
                         (True if edge.is_internal() and edge.head_node
                         != seed and edge.tail_node
                         != seed else False)))

        # Test whether root edge is eligible (i.e., the edge is not a
        # leaf when the tree is unrooted). If the test passes, add 'root'
        # tp the list of eligible edges

        if not any([x.is_leaf() for x in seed.child_nodes()]):
            edge_list += ['root']

        chosen_edge = random.choice(edge_list)
        print chosen_edge

        # The NNI is done with the chosen edge as root. This will
        # involve rerooting the tree if we choose an edge that isn't
        # the root

        if chosen_edge != 'root':
            tree.reroot_at_edge(chosen_edge, length1=chosen_edge.length
                                / 2, length2=chosen_edge.length / 2,
                                delete_outdegree_one=False)
            root = tree.seed_node
        else:
            root = seed

        # To do the swap: find the nodes on either side of root

        (child_left, child_right) = root.child_nodes()

        # Pick a child node from each of these

        neighbour1 = random.choice(child_left.child_nodes())
        neighbour2 = random.choice(child_right.child_nodes())

        # Prune the chosen nearest neighbours - but don't
        # do any tree structure updates

        tree.prune_subtree(neighbour1, update_splits=False,
                           delete_outdegree_one=False)
        tree.prune_subtree(neighbour2, update_splits=False,
                           delete_outdegree_one=False)

        # Reattach the pruned neighbours to the opposite side
        # of the tree

        child_left.add_child(neighbour2)
        child_right.add_child(neighbour1)

        # Reroot the tree using the original seed node, and
        # update splits

        if not chosen_edge == 'root':
            tree.reroot_at_node(seed, update_splits=True)
        else:
            tree.update_splits()

        if resolved:
            print 'Reinstating root trisomy'
            tree.deroot()

        newick = tree.as_newick_string()
        if not newick.endswith(';'):
            newick += ';'
        if tree.is_rooted:
            newick = '[&R] ' + newick

        return Tree(newick, self.score, self.program, self.name,
                    self.output)

    def spr(self, disallow_sibling_SPRs=False):

        def choose_prune_and_regraft_nodes(time, blocks,
                disallow_sibling_SPRs):
            matching_branches = [x for x in dists if x[1] < time
                                 < x[2]]

            prune = random.sample(matching_branches, 1)[0]

            if disallow_sibling_SPRs:
                siblings = prune[0].sister_nodes()
                for br in matching_branches:
                    if br[0] in siblings:
                        matching_branches.remove(br)

            matching_branches.remove(prune)
            
            if matching_branches == []:
                print 'No non-sibling branches available'
                return None, None

            regraft = random.sample(matching_branches, 1)[0]

            prune_taxa = [n.taxon.label for n in prune[0].leaf_iter()]
            regraft_taxa = [n.taxon.label for n in regraft[0].leaf_iter()]
            print 'Donor group = {0}'.format(regraft_taxa)
            print 'Receiver group = {0}'.format(prune_taxa)
            return (prune, regraft)

        def get_blocks(tree, include_leaf_nodes=True):
            dists = []
            blocks = {}
            if include_leaf_nodes:
                iterator = tree.preorder_node_iter()
            else:
                iterator = tree.preorder_internal_node_iter()
            for n in iterator:
                node_height = n.distance_from_root()
                if not n.parent_node:
                    root_height = n.distance_from_root()
                    tree_height = root_height + n.distance_from_tip()
                    parent_height = 0
                else:
                    parent_height = n.parent_node.distance_from_root()
                node_height = round(node_height, 8)
                parent_height = round(parent_height, 8)

                if not node_height in blocks:
                    blocks[node_height] = []

                dists.append((n, parent_height, node_height))

            for time in blocks:
                for (node, parent_h, node_h) in dists:
                    if parent_h < time <= node_h:
                        blocks[time].append(node)

            dists.sort(key=lambda x: x[2])
            return (blocks, dists)

        def weight_by_branches(blocks):
            intervals = sorted(blocks.keys())
            weighted_intervals = [0] + [None] * (len(intervals) - 1)
            for i in range(1, len(intervals)):
                time_range = intervals[i] - intervals[i - 1]
                num_branches = len(blocks[intervals[i]])
                weighted_range = time_range * num_branches
                weighted_intervals[i] = weighted_range \
                    + weighted_intervals[i - 1]
            return weighted_intervals

        def get_time(blocks, weights=None):
            d = sorted(blocks.keys())
            if weights:
                samp = random.uniform(weights[0], weights[-1])
                for i in range(len(weights) - 1):
                    if weights[i + 1] >= samp > weights[i]:
                        interval = weights[i + 1] - weights[i]
                        proportion = (samp - weights[i]) / interval
                        break
                drange = d[i + 1] - d[i]
                time = drange * proportion + d[i]
            else:
                time = random.uniform(d[0], d[-1])

            print 'LGT event at time: {0}'.format(time)

            return time

        def add_node(tree, time, regraft_node):
            parent_node = regraft_node[0].parent_node
            new_node = parent_node.add_child(dpy.Node(),
                    edge_length=time - regraft_node[1])
            tree.reindex_subcomponent_taxa()
            tree.update_splits()
            return new_node

        def prunef(tree, node):
            tree.prune_subtree(node, update_splits=False,
                               delete_outdegree_one=True)

        def regraftf(
            tree,
            time,
            target_node,
            child_node,
            ):

            target_node.add_child(child_node[0],
                                  edge_length=child_node[2] - time)
            return tree

        tree = dpy.Tree()
        tree.read_from_string(self.newick, 'newick')
        (blocks, dists) = get_blocks(tree)
        time = get_time(blocks)
        (p, r) = choose_prune_and_regraft_nodes(time, dists,
                disallow_sibling_SPRs=disallow_sibling_SPRs)
        if (p, r) == (None, None):
            return self
        new_node = add_node(tree, time, r)
        prunef(tree, p[0])
        prunef(tree, r[0])
        regraftf(tree, time, new_node, p)
        regraftf(tree, time, new_node, r)
        tree.reindex_subcomponent_taxa()
        tree.update_splits()

        newick = tree.as_newick_string()
        if not newick.endswith(';'):
            newick += ';'

        if tree.is_rooted:
            newick = '[&R] ' + newick

        return Tree(newick, self.score, self.program, self.name,
                    self.output)
