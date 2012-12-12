#!/usr/bin/env python

from tree import Tree
from sequence_record import TCSeqRec
from errors import directorycheck_and_make, directorycheck_and_quit

# import GeoMeTreeHack

import dendropy as dpy
import glob
import numpy as np
import os
import random
import re
import shutil
import sys
import copy
from textwrap import dedent

np.set_printoptions(precision=3, linewidth=200)


class SeqSim(object):

    """ 
    This class is a front end to the ALF simulator
    """

    def __init__(
        self,
        simulation_name='sim',
        working_directory='./alftmp/',
        outfile_path='./',
        unit_is_pam=True,
        ):

        if unit_is_pam is True:
            unit_is_pam = 'true'
        else:
            unit_is_pam = 'false'

        file_string = \
            dedent('''\
            # Filepath parameters
            ##############################
            # directories for file storage
            wdir := '{0}';
            dbdir := 'DB/';
            dbAncdir := 'DBancestral/';

            '''.format(working_directory))

        pam_string = \
            dedent('''\
            # Time units
            ##########################
            # timescale for simulation
            unitIsPam := {0}:

            '''.format(unit_is_pam))

        name_string = \
            dedent('''\
            # Simulation Name
            ####################
            # name of simulation
            mname := {0};

            '''.format(simulation_name))

        self.outfile_path = outfile_path
        self.name = simulation_name

        self.parameters = {
            'files': file_string,
            'pam': pam_string,
            'name': name_string,
            'genome': '',
            'subst': '',
            'indels': '',
            'ratevar': '',
            'tree': '',
            }

    def __str__(self):

        sections = []
        for section in [
            'files',
            'pam',
            'name',
            'genome',
            'subst',
            'indels',
            'ratevar',
            'tree',
            ]:
            if section in self.parameters:
                sections.append(self.parameters[section])
        return ''.join(sections)

    def filepaths(self, simulation_name='sim', working_directory='./'):

        file_string = \
            dedent('''\
            # Filepath parameters
            ##############################
            # directories for file storage
            wdir := '{0}';
            dbdir := 'DB/';
            dbAncdir := 'DBancestral/';

            # name of simulation
            mname := {1};

            '''.format(working_directory,
                   simulation_name))

        self.parameters['files'] = file_string
        return file_string

    def root_genome(
        self,
        number_of_genes=100,
        min_length=10,
        kappa=1,
        theta=1,
        ):

        genome_string = \
            dedent('''\
            # Root genome parameters
            #############################
            realseed := false;
            protStart := {0};
            minGeneLength := {1};
            gammaLengthDist := [{2},{3}];
            blocksize := 1:

            '''.format(number_of_genes,
                   min_length, kappa, theta))

        self.parameters['genome'] = genome_string
        return genome_string

    def pam(self, unit_is_pam=True):

        if unit_is_pam is True:
            unit_is_pam = 'true'
        else:
            unit_is_pam = 'false'

        pam_string = \
            dedent('''\
            # Time units
            ##########################
            # timescale for simulation
            unitIsPam := {0}:

            '''.format(unit_is_pam))

        self.parameters['pam'] = pam_string
        return pam_string

    def rename(self, name):

        name_string = \
            dedent('''\
            # Simulation Name
            ####################
            # name of simulation
            mname := {0};

            '''.format(name))

        self.parameters['name'] = name_string
        self.name = name
        return name_string

    def gtr_model(
        self,
        CtoT,
        AtoT,
        GtoT,
        AtoC,
        CtoG,
        AtoG,
        Afreq,
        Cfreq,
        Gfreq,
        Tfreq,
        allow_nonsense=False,
        ):

        if allow_nonsense:
            allow_nonsense = 'true'
        else:
            allow_nonsense = 'false'

        subst_string = \
            dedent('''\
            # Substitution Model Parameters
            ######################################################################################################
            substModels := [SubstitutionModel('GTR', [{0}, {1}, {2}, {3}, {4}, {5}], [{6}, {7}, {8}, {9}], {10})];

            '''.format(
            CtoT,
            AtoT,
            GtoT,
            AtoC,
            CtoG,
            AtoG,
            Afreq,
            Cfreq,
            Gfreq,
            Tfreq,
            allow_nonsense,
            ))

        self.parameters['subst'] = subst_string
        return subst_string

    def hky_model(
        self,
        alpha,
        beta,
        Afreq,
        Cfreq,
        Gfreq,
        Tfreq,
        allow_nonsense=False,
        ):

        if allow_nonsense:
            allow_nonsense = 'true'
        else:
            allow_nonsense = 'false'

        subst_string = \
            dedent('''\
            # Substitution Model Parameters
            #################################################################################
            substModels := [SubstitutionModel('HKY', [{0}, {1}], [{2}, {3}, {4}, {5}], {6})];

            '''.format(
            alpha,
            beta,
            Afreq,
            Cfreq,
            Gfreq,
            Tfreq,
            allow_nonsense,
            ))

        self.parameters['subst'] = subst_string
        return subst_string

    def f84_model(
        self,
        kappa,
        beta,
        Afreq,
        Cfreq,
        Gfreq,
        Tfreq,
        allow_nonsense=False,
        ):

        if allow_nonsense:
            allow_nonsense = 'true'
        else:
            allow_nonsense = 'false'

        subst_string = \
            dedent('''\
            # Substitution Model Parameters
            #################################################################################
            substModels := [SubstitutionModel('F84', [{0}, {1}], [{2}, {3}, {4}, {5}], {6})];

            '''.format(
            kappa,
            beta,
            Afreq,
            Cfreq,
            Gfreq,
            Tfreq,
            allow_nonsense,
            ))
        return subst_string

    def jc_model(self, allow_nonsense=False):

        if allow_nonsense:
            allow_nonsense = 'true'
        else:
            allow_nonsense = 'false'

        subst_string = \
            dedent('''\
            # Substitution Model Parameters
            ######################################################################
            substModels := [SubstitutionModel('F84', [1, 1], [seq(0.25,4)], {0})];

            '''.format(allow_nonsense))

        self.parameters['subst'] = subst_string
        return subst_string

    def one_word_model(self, model='WAG'):
        """ 
        "One-word" models don't need any extra parameters.
        Codon: CPAM, ECM, ECMu
        AA:    GCB,  JTT, LG,  WAG
        """

        subst_string = \
            dedent('''\
            # Substitution Model Parameters
            ##########################################
            substModels := [SubstitutionModel('{0}')];

            '''.format(model))

        self.parameters['subst'] = subst_string
        return subst_string

    def indels(
        self,
        gain_rate=0.0001,
        gain_model='ZIPF',
        gain_params=[1.821],
        max_gain_length=10,
        loss_rate=0.0001,
        loss_model='ZIPF',
        loss_params=[1.821],
        max_loss_length=10,
        ):

        indel_string = \
            dedent('''\
            # Indel Parameters
            ####################################################################
            indelModels := [IndelModel({0}, {1}, {2}, {3}, {4}, {5}, {6}, {7})];

            '''.format(
            gain_rate,
            gain_model,
            gain_params,
            max_gain_length,
            loss_rate,
            loss_model,
            loss_params,
            max_loss_length,
            ))

        self.parameters['indels'] = indel_string
        return indel_string

    def rate_variation(
        self,
        shape=1,
        ncat=4,
        pinvar=0,
        ):
        """
        Models rate variation among sites
        Only uses gamma distribution (ALF allows more options)
        """

        rate_var_string = \
            dedent('''\
            # Rate Variation Parameters
            ########################################################
            rateVarModels := [RateVarModel(Gamma, {0}, {1}, {2})];

            '''.format(ncat,
                   pinvar, shape))

        self.parameters['ratevar'] = rate_var_string
        return rate_var_string

    def custom_tree(self, treefile):

        custom_tree_string = \
            dedent('''\
            # Tree Parameters
            #####################
            treeType := 'Custom';
            treeFile := '{0}';

            '''.format(treefile))

        self.parameters['tree'] = custom_tree_string
        return custom_tree_string

    def BDtree(
        self,
        birthrate=0.01,
        deathrate=0.001,
        nspecies=15,
        mutation_rate=250,
        scale_tree=True,
        ultrametric=False,
        ):

        if ultrametric:
            ultrametric = 'true'
        else:
            ultrametric = 'false'

        if scale_tree:
            scale_tree = 'true'
        else:
            scale_tree = 'false'

        bdtree_string = \
            dedent('''\
            # Tree Parameters
            #####################
            treeType := 'BDTree';
            birthRate := {0};
            deathRate := {1};
            NSpecies := {2};
            ultrametric := {3};
            mutRate := {4};
            scaleTree := {5};

            '''.format(
            birthrate,
            deathrate,
            nspecies,
            ultrametric,
            mutation_rate,
            scale_tree,
            ))

        self.parameters['tree'] = bdtree_string
        return bdtree_string

    def write_parameters(self):
        outfile_name = '{0}/{1}.drw'.format(self.outfile_path,
                self.name)
        writer = open(outfile_name, 'w')
        writer.write(self.__str__())
        writer.flush()
        writer.close()

        return outfile_name

    def runALF(self, parameter_file, quiet=True):
        """
        alfsim must be in path
        """

        if not os.path.isfile(parameter_file):
            self.write_parameters(parameter_file)

        print 'Running ALF on {0}'.format(parameter_file)
        command = 'alfsim {0}'.format(parameter_file)
        if quiet:
            command += ' > /dev/null 2> /dev/null'
        os.system(command)

        return

    def check_diff_top(self, check_tree, tree_list):
        """
        Returns True if topologies are different,
        False if they are the same (unweighted RF
        distance = 0)
        """

        checklist = []
        t1 = dpy.Tree()
        t1.read_from_string(check_tree.newick, 'newick')

        for tree in tree_list:
            t2 = dpy.Tree()
            t2.read_from_string(tree.newick, 'newick')
            if t2.symmetric_difference(t1):
                checklist.append(True)
            else:
                checklist.append(False)

        if all(checklist):
            return True  # evaluates to True if list is empty
        else:

            return False

    @classmethod
    def simulate_from_tree_GTR(
        cls,
        tree,
        name,
        tmpdir='/tmp',
        allow_nonsense=False,
        ):
        """
        parameters = dictionary with keys:
                     Afreq, Cfreq, Gfreq, Tfreq, 
                     AtoC, AtoG, AtoT, CtoG, CtoT, GtoT,
                     gamma.
        tree       = Tree object
        """

        directorycheck_and_quit(tmpdir)
        GTR_parameters = tree.extract_GTR_parameters()
        gamma = tree.extract_gamma_parameter()
        param_dir = '{0}/alf_parameter_dir'.format(tmpdir)
        working_dir = '{0}/alf_working_dir'.format(tmpdir)
        treefile = '{0}/treefile.nwk'.format(tmpdir)

        tree.pam2sps('sps2pam').write_to_file(treefile)

        directorycheck_and_make(param_dir)
        directorycheck_and_make(working_dir)

        sim = cls(simulation_name=name, working_directory=working_dir,
                  outfile_path=param_dir, unit_is_pam=True)

        sim.indels()
        sim.rate_variation(gamma)
        sim.root_genome(number_of_genes=1, min_length=length)
        sim.gtr_model(
            CtoT=GTR_parameters['CtoT'],
            AtoT=GTR_parameters['AtoT'],
            GtoT=GTR_parameters['GtoT'],
            AtoC=GTR_parameters['AtoC'],
            CtoG=GTR_parameters['CtoG'],
            AtoG=GTR_parameters['AtoG'],
            Afreq=GTR_parameters['Afreq'],
            Cfreq=GTR_parameters['Cfreq'],
            Gfreq=GTR_parameters['Gfreq'],
            Tfreq=GTR_parameters['Tfreq'],
            )
        sim.custom_tree(treefile)
        params = sim.write_parameters()
        sim.runALF(params, quiet=True)
        tree_newick = tree.newick
        alf_newick = \
            open('{0}/alf_working_dir/{1}/RealTree.nwk'.format(tmpdir,
                 name)).read()
        replacement_dict = dict(zip(re.findall(r'(\w+)(?=:)',
                                alf_newick), re.findall(r'(\w+)(?=:)',
                                tree_newick)))  # bug correction

    @classmethod
    def simulate_from_tree_WAG(
        cls,
        self,
        tree,
        gamma,
        simulation_name='sim',
        working_directory='./alftmp/',
        outfile_path='./',
        unit_is_pam=True,
        tmpdir='/tmp',
        ):

        sim = cls(simulation_name, working_directory, outfile_path,
                  unit_is_pam)
        sim.indels()
        sim.rate_variation(gamma)
        sim.root_genome(number_of_genes=1, min_length=length)
        sim.one_word_model('WAG')
        sim.custom_tree(treefile)
        sim.write_parameters()
        sim.runALF()
        tree_newick = open(treefile).read()
        alf_newick = \
            open('{0}/alf_working_dir/{1}/RealTree.nwk'.format(tmpdir,
                 name)).read()
        replacement_dict = dict(zip(re.findall(r'(\w+)(?=:)',
                                alf_newick), re.findall(r'(\w+)(?=:)',
                                tree_newick)))  # bug correction

    def simulate_set(
        self,
        K,
        M,
        n,
        tune,
        regime,
        branch_length_func,
        inner_edge_params,
        leaf_params,
        scale_func,
        master_tree_generator_method='random_topology',
        class_tree_permuter='nni',
        guarantee_unique=False,
        num_permutations=0,
        scale_params=(1, 1),
        gene_length_kappa=1,
        gene_length_theta=1,
        gene_length_min=10,
        filepath='./',
        tmpdir='/tmp',
        gtp_path='./class_files',
        unit_is_pam=True,
        quiet=True,
        ):
        """
        Regime 1:
            1 topology (n species)
            M alignments
            (2n - 3) branch lengths in total

        Regime 2:
            K topologies (n species)
            M alignments, distributed among K classes
            K * (2n - 3) branch lengths in total

        Regime 3:
            K topologies (n species)
            M alignments, distributed among K classes
            Each of Mk alignments in class k has scaled branch lengths 
            (Mk - 1) * (2n - 3) branch lengths in total

        Regime 4:
            K topologies (n species)
            M alignments, distributed among K classes
            Each of Mk alignments in class k has independent branch lengths
            M * K * (2n - 3) branch lengths in total

        Tuning:
        The tuning parameter gives coarse control over the difference in sizes
        of groups - for example a very large value ( > 1000 ) tends to
        assign M - K + 1 genes to a single group, and 1 gene to each of the
        others, and a very small value ( < 1/1000 ) tends to assign M/K genes
        to each class. A zero value makes all groups the same size.

        """

        print '{0} = {1}'.format(class_tree_permuter, num_permutations)

        def class_stats(M, mk_list):
            d = {}
            nclasses = len(mk_list)
            Msize = len(M)
            ind = np.triu_indices(Msize, 1)
            intra_class = []
            inter_class = []
            cs = np.concatenate((np.array([0]), np.cumsum(mk_list)))
            for i in range(nclasses):
                intra_class += list(M[cs[i]:cs[i + 1], cs[i]:cs[i
                                    + 1]][np.triu_indices(mk_list[i],
                                    1)].flatten())
                inter_class += list(M[cs[i]:cs[i + 1], cs[i + 1]:
                                    ].flatten())
            d['overall_mean'] = np.mean(M[ind])
            d['intra_mean'] = np.mean(intra_class)
            d['inter_mean'] = np.mean(inter_class)
            d['overall_var'] = np.var(M[ind])
            d['intra_var'] = np.var(intra_class)
            d['inter_var'] = np.var(inter_class)

            return d

        def make_master_tree(
            n,
            names,
            method,
            inner_edge_params=(1, 1),
            leaf_params=(1, 1),
            distribution_func=np.random.gamma,
            ):
            """
            Function returns a tree object with n tips,
            named according to `names`, and constructed
            according to `method`, which is one of 'random_topology',
            'random_yule' and 'random_coal'
            """

            if method == 'random_topology':
                master_topology = Tree().random_topology(n,
                        names=names, rooted=True)
                master_tree = \
                    master_topology.randomise_branch_lengths(inner_edges=inner_edge_params,
                        leaves=leaf_params,
                        distribution_func=branch_length_func)
                master_tree.newick = '[&R] ' + master_tree.newick
            elif method == 'random_yule':
                master_tree = Tree().random_yule(n, names=names)
            elif method == 'random_coal':
                master_tree = Tree().random_coal(n, names=names)
            return master_tree

        def make_class_tree(
            master_tree,
            permutation_extent,
            method,
            with_check=True,
            checklist=[],
            ):
            """
            Function returns a tree object derived from a master tree,
            but with some permutation applied. The type of permutation
            is defined by `method`, one of 'nni', 'spr' and 'coal'
            If with_check is True, the generated tree is checked against
            a checklist of other trees on the same species, and permutations
            are applied until the new tree has a unique topology. This is 
            only implemented for nni and spr.
            """

            if num_permutations == 0:
                return master_tree

            new_tree = Tree(master_tree.newick)

            if method == 'nni':
                if with_check:
                    while not self.check_diff_top(new_tree, checklist):
                        new_tree = Tree(master_tree.newick)
                        for i in range(permutation_extent):
                            new_tree = new_tree.nni()
                else:
                    for i in range(num_permutations):
                        new_tree = new_tree.nni()
            elif method == 'spr':
                if with_check:
                    while not self.check_diff_top(new_tree, checklist):
                        new_tree = Tree(master_tree.newick)
                        for i in range(permutation_extent):
                            new_tree = \
                                new_tree.spr(disallow_sibling_SPRS=True)
                else:
                    for i in range(num_permutations):
                        new_tree = new_tree.spr()
            elif method == 'coal':
                new_tree = \
                    master_tree.get_constrained_gene_tree(scale_to=permutation_extent)
            return new_tree

        # Create directories for simulation trees and parameter files

        if not os.path.isdir('{0}/alf_parameter_dir'.format(tmpdir)):
            os.mkdir('{0}/alf_parameter_dir'.format(tmpdir))
        if not os.path.isdir('{0}/alf_trees_dir'.format(tmpdir)):
            os.mkdir('{0}/alf_trees_dir'.format(tmpdir))
        if not os.path.isdir(filepath):
            os.mkdir(filepath)
        if not os.path.isdir('{0}/true_trees'.format(filepath)):
            os.mkdir('{0}/true_trees'.format(filepath))
        if not os.path.isdir('{0}/true_trees/individual'.format(filepath)):
            os.mkdir('{0}/true_trees/individual'.format(filepath))
        if not os.path.isdir('{0}/dna_alignments'.format(filepath)):
            os.mkdir('{0}/dna_alignments'.format(filepath))
        if not os.path.isdir('{0}/aa_alignments'.format(filepath)):
            os.mkdir('{0}/aa_alignments'.format(filepath))

        # Assign numbers of genes to classes
        # list `mk` gives number of genes in each class

        if regime == 1:
            K = 1

        if tune == 0:
            proportions = [float(K) / M for x in range(K)]
        else:

            proportions = np.random.gamma(shape=float(M) / (tune * K),
                    scale=tune * float(K) / M, size=K)

        s = sum(proportions)
        mk = [int((np.round(x * M / s) if x * M / s > 0.5 else 1.0))
              for x in proportions]
        diff = M - sum(mk)
        if diff > 0:
            mk[mk.index(min(mk))] += diff
        else:
            mk[mk.index(max(mk))] += diff
        assert min(mk) > 0.0

        true_clustering = []
        for i in range(K):
            for j in range(mk[i]):
                true_clustering.append(i + 1)

        print 'Simulating {0} genes in {1} classes, distributed as {2}'.format(M,
                K, mk)

        names = ['Sp{0}'.format(i) for i in range(1, n + 1)]

        print 'N classes =', K
        print 'N genes = ', M
        print 'N species =', n
        print 'Regime = ', regime
        print 'N permutations =', num_permutations
        print 'Tuning =', tune
        print 'mk =', mk
        print 'true clustering = ', true_clustering

        # Create simulation trees

        # Make a master tree

        master_tree = make_master_tree(
            n,
            names,
            method=master_tree_generator_method,
            inner_edge_params=inner_edge_params,
            leaf_params=leaf_params,
            distribution_func=branch_length_func,
            )
        class_trees = []
        parameter_files = []

        print 'Master tree = ', master_tree
        master_tree.write_to_file('{0}/true_trees/master.tree'.format(filepath),
                                  suppress_NHX=True)

        # make K class trees

        for k in range(K):
            print 'Making class {0}/{1}'.format(k + 1, K)

            if num_permutations > 0:
                class_tree = make_class_tree(master_tree,
                        num_permutations, class_tree_permuter,
                        with_check=guarantee_unique,
                        checklist=class_trees)
                class_trees.append(class_tree)
            else:

                class_tree = make_master_tree(
                    n,
                    names,
                    method=master_tree_generator_method,
                    inner_edge_params=inner_edge_params,
                    leaf_params=leaf_params,
                    distribution_func=branch_length_func,
                    )
                class_trees.append(class_tree)

            print 'class tree = ', class_tree
            class_tree.write_to_file('{0}/true_trees/class{1}.tree'.format(filepath,
                    k + 1), suppress_NHX=True)

            # ALF only behaves itself if trees are in PAM units,
            # so we scale our newly-generated class trees to have branch lengths
            # in PAM units.
            # Our class_trees list contains unconverted trees

            if unit_is_pam:  # Default = True
                class_tree_PAM = class_tree.pam2sps('sps2pam')  # conversion from SPS to PAM
                class_tree_PAM.write_to_file('{0}/alf_trees_dir/class{1}_1.nwk'.format(tmpdir,
                        k + 1), suppress_NHX=True)

            # Write parameter files

            ngenes = mk[k]

            sim = SeqSim(simulation_name='class{0}_1'.format(k + 1),
                         working_directory='{0}/alf_working_dir'.format(tmpdir),
                         outfile_path='{0}/alf_parameter_dir'.format(tmpdir),
                         unit_is_pam=unit_is_pam)  # make new simulation object
            sim.parameters['subst'] = self.parameters['subst']  # copy over global parameters
            sim.parameters['indels'] = self.parameters['indels']
            sim.parameters['ratevar'] = self.parameters['ratevar']

            if regime in [1, 2]:

                sim.root_genome(number_of_genes=ngenes,
                                kappa=gene_length_kappa,
                                theta=gene_length_theta)
                sim.custom_tree('{0}/alf_trees_dir/class{1}_1.nwk'.format(tmpdir,
                                k + 1))
                sim.write_parameters()

                continue

            # For regimes 3 & 4 each gene within a class is simulated along its own tree:
            #   Under regime 3 each gene within a class has its branch lengths scaled
            #   by a scaling parameter, and under regime 4 each gene has random branch lengths
            #   drawn from some distribution (although individual branch rescaling might be better)

            for genes in range(ngenes):
                if regime == 3:
                    scale_factor = scale_func(*scale_params)
                    class_tree = class_trees[k].scale(scale_factor)
                elif regime == 4:
                    class_tree = \
                        class_trees[k].randomise_branch_lengths(inner_edges=inner_edge_params,
                            leaves=leaf_params,
                            distribution_func=branch_length_func)
                    if unit_is_pam:  # same conversion as before
                        class_tree = class_tree.pam2sps('sps2pam')

                class_tree.write_to_file('{0}/alf_trees_dir/class{1}_{2}.nwk'.format(tmpdir,
                        k + 1, genes + 1), suppress_NHX=True)

                sim.root_genome(number_of_genes=1,
                                min_length=gene_length_min,
                                kappa=gene_length_kappa,
                                theta=gene_length_theta)
                sim.custom_tree('{0}/alf_trees_dir/class{1}_{2}.nwk'.format(tmpdir,
                                k + 1, genes + 1))
                sim.rename('class{0}_{1}'.format(k + 1, genes + 1))
                sim.write_parameters()

        # Estimate distances between base class trees

        # if unit_is_pam:
        #     class_trees = [x.pam2sps() for x in class_trees]

        geodists = []
        eucdists = []
        symdists = []
        wrfdists = []
        with open('{0}/basetrees.nwk'.format(tmpdir), 'w') as file:
            file.write('\n'.join([x.newick.rstrip() for x in
                       class_trees]))
        os.system('java -jar {0}/gtp.jar -o {1}/baseout.txt {1}/basetrees.nwk'.format(gtp_path,
                  tmpdir))
        with open('{0}/baseout.txt'.format(tmpdir)) as file:
            for line in file:
                line = line.rstrip()
                if line:
                    (i, j, value) = line.split()
                    geodists.append(float(value))
        for a in range(K):
            tree_a = dpy.Tree.get_from_string(class_trees[a].newick,
                    'newick')
            for b in range(a + 1, K):
                tree_b = \
                    dpy.Tree.get_from_string(class_trees[b].newick,
                        'newick')
                eucdists.append(tree_a.euclidean_distance(tree_b))
                symdists.append(tree_a.symmetric_difference(tree_b))
                wrfdists.append(tree_a.robinson_foulds_distance(tree_b))

        writer = open('{0}/treedistances.txt'.format(filepath), 'w')
        writer.write('''True clustering:\t{0}
Class base tree distances:
geodesic\t{1}
euclidean\t{2}
RF\t{3}
wRF\t{4}

'''.format(true_clustering,
                     np.mean(geodists), np.mean(eucdists),
                     np.mean(symdists), np.mean(wrfdists)))
        writer.flush()

        # Run simulations, and correct ALF renaming bug

        parameter_files = \
            glob.glob('{0}/alf_parameter_dir/*.drw'.format(tmpdir))
        tree_files = glob.glob('{0}/alf_trees_dir/*.nwk'.format(tmpdir))
        sort_key = lambda item: tuple((int(num) if num else alpha)
                for (num, alpha) in re.findall(r'(\d+)|(\D+)', item))
        parameter_files.sort(key=sort_key)
        tree_files.sort(key=sort_key)
        files = zip(parameter_files, tree_files)

        for (params, tree) in files:
            self.runALF(params, quiet=quiet)
            name = params[params.rindex('/'):params.rindex('.')]
            (class_number, base_gene_number) = re.findall(r'(\d+)',
                    name)
            tree_newick = open(tree).read()
            alf_newick = \
                open('{0}/alf_working_dir/{1}/RealTree.nwk'.format(tmpdir,
                     name)).read()
            replacement_dict = dict(zip(re.findall(r'(\w+)(?=:)',
                                    alf_newick),
                                    re.findall(r'(\w+)(?=:)',
                                    tree_newick)))  # bug correction

            for dna_alignment in \
                sorted(glob.glob('{0}/alf_working_dir/{1}/MSA/*dna.fa'.format(tmpdir,
                       name)), key=sort_key):
                gene_number = dna_alignment[dna_alignment.rindex('/')
                    + 1:].split('_')[1]
                record = TCSeqRec(dna_alignment)
                record.sort_by_name()
                record.headers = [replacement_dict[x[:x.rindex('/')]]
                                  for x in record.headers]
                record.write_fasta('{0}/dna_alignments/class{1}_{2}.fas'.format(filepath,
                                   class_number, int(base_gene_number)
                                   + int(gene_number) - 1))
                record.write_phylip('{0}/dna_alignments/class{1}_{2}.phy'.format(filepath,
                                    class_number, int(base_gene_number)
                                    + int(gene_number) - 1))

            for aa_alignment in \
                sorted(glob.glob('{0}/alf_working_dir/{1}/MSA/*aa.fa'.format(tmpdir,
                       name)), key=sort_key):
                gene_number = aa_alignment[aa_alignment.rindex('/')
                    + 1:].split('_')[1]
                record = TCSeqRec(aa_alignment)
                record.sort_by_name()
                record.headers = [replacement_dict[x[:x.rindex('/')]]
                                  for x in record.headers]
                record.write_fasta('{0}/aa_alignments/class{1}_{2}.fas'.format(filepath,
                                   class_number, int(base_gene_number)
                                   + int(gene_number) - 1))
                record.write_phylip('{0}/aa_alignments/class{1}_{2}.phy'.format(filepath,
                                    class_number, int(base_gene_number)
                                    + int(gene_number) - 1))

            # Write true trees

            if regime in [1, 2]:
                for g in range(mk[int(class_number) - 1]):
                    Tree(tree_newick).pam2sps().write_to_file('{0}/true_trees/individual/class{1}_{2}.nwk'.format(filepath,
                            class_number, g + 1), suppress_NHX=True)
            else:

                Tree(tree_newick).pam2sps().write_to_file('{0}/true_trees/individual/{1}.nwk'.format(filepath,
                        name), suppress_NHX=True)

        # Intra- and inter-class stats

        alltrees = \
            glob.glob('{0}/true_trees/individual/*.nwk'.format(filepath))

        alltrees.sort(key=sort_key)

        alltrees = [open(x).read().rstrip() for x in alltrees]

        dpytrees = [dpy.Tree.get_from_string(x, 'newick') for x in
                    alltrees]

        # for x in range(len(alltrees)):
        #     print x,'\n',alltrees[x], '\n',dpy.Tree.get_from_string(alltrees[x],'newick').as_newick_string()

        geodists = np.zeros([M, M])
        eucdists = np.zeros([M, M])
        symdists = np.zeros([M, M])
        wrfdists = np.zeros([M, M])

        # using gtp.jar for geodesic distances

        with open('{0}/geotrees.nwk'.format(tmpdir), 'w') as file:
            file.write('\n'.join(alltrees))
        os.system('java -jar {0}/gtp.jar -o {1}/output.txt {1}/geotrees.nwk'.format(gtp_path,
                  tmpdir))
        with open('{0}/output.txt'.format(tmpdir)) as file:
            for line in file:
                line = line.rstrip()
                if line:
                    (i, j, value) = line.split()
                    i = int(i)
                    j = int(j)
                    value = float(value)
                    geodists[i, j] = geodists[j, i] = value

        for a in range(M):
            for b in range(a + 1, M):
                eucdists[a, b] = eucdists[b, a] = \
                    dpytrees[a].euclidean_distance(dpytrees[b])
                symdists[a, b] = symdists[b, a] = \
                    dpytrees[a].symmetric_difference(dpytrees[b])
                wrfdists[a, b] = wrfdists[b, a] = \
                    dpytrees[a].robinson_foulds_distance(dpytrees[b])

        geodic = class_stats(geodists, mk)
        eucdic = class_stats(eucdists, mk)
        symdic = class_stats(symdists, mk)
        wrfdic = class_stats(wrfdists, mk)

        writer.write('Geodesic class stats\n')
        for key in sorted(geodic):
            writer.write('{0}\t{1}\n'.format(key, geodic[key]))
        writer.write('\n')
        writer.flush()

        writer.write('Euc class stats\n')
        for key in sorted(eucdic):
            writer.write('{0}\t{1}\n'.format(key, eucdic[key]))
        writer.write('\n')
        writer.flush()

        writer.write('RF class stats\n')
        for key in sorted(symdic):
            writer.write('{0}\t{1}\n'.format(key, symdic[key]))
        writer.write('\n')
        writer.flush()

        writer.write('wRF class stats\n')
        for key in sorted(wrfdic):
            writer.write('{0}\t{1}\n'.format(key, wrfdic[key]))
        writer.write('\n')
        writer.flush()

        writer.close()

        shutil.rmtree('{0}/alf_parameter_dir'.format(tmpdir))
        shutil.rmtree('{0}/alf_trees_dir'.format(tmpdir))
        shutil.rmtree('{0}/alf_working_dir'.format(tmpdir))
        os.remove('{0}/output.txt'.format(tmpdir))
        os.remove('{0}/geotrees.nwk'.format(tmpdir))
        os.remove('{0}/basetrees.nwk'.format(tmpdir))
        os.remove('{0}/baseout.txt'.format(tmpdir))


if __name__ == '__main__':

    # Small test suite

    sim = SeqSim()
    sim.hky_model(
        alpha=3.551,
        beta=1,
        Afreq=0.31,
        Cfreq=0.18,
        Gfreq=0.21,
        Tfreq=0.3,
        )
    sim.indels()
    sim.rate_variation()

    # sim.simulate_set(K=4,M=20,n=12,tune=1,regime=2,branch_length_func=np.random.gamma,inner_edge_params=(3.2,0.029),leaf_params=(2.2,0.097),scale_func=np.random.gamma,
    #     master_tree_generator_method='random_yule',class_tree_permuter='spr',guarantee_unique=True,num_permutations=3,scale_params=(1,1),gene_length_kappa=5.53,gene_length_theta=72.35,gene_length_min=10,
    #     filepath='/Users/kgori/scratch/testspr',tmpdir='/tmp',gtp_path='/Users/kgori/git/kevin/clustering_project/class_files/',unit_is_pam=True,quiet=True)

    # sim.simulate_set(K=4,M=20,n=12,tune=1,regime=2,branch_length_func=np.random.gamma,inner_edge_params=(3.2,0.029),leaf_params=(2.2,0.097),scale_func=np.random.gamma,
    #     master_tree_generator_method='random_yule',class_tree_permuter='nni',guarantee_unique=True,num_permutations=1,scale_params=(1,1),gene_length_kappa=5.53,gene_length_theta=72.35,gene_length_min=10,
    #     filepath='/Users/kgori/scratch/testnni',tmpdir='/tmp',gtp_path='/Users/kgori/git/kevin/clustering_project/class_files/',unit_is_pam=True,quiet=True)

    # sim.simulate_set(K=4,M=20,n=12,tune=1,regime=2,branch_length_func=np.random.gamma,inner_edge_params=(3.2,0.029),leaf_params=(2.2,0.097),scale_func=np.random.gamma,
    #     master_tree_generator_method='random_yule',class_tree_permuter='coal',guarantee_unique=True,num_permutations=200,scale_params=(1,1),gene_length_kappa=5,gene_length_theta=100,gene_length_min=25,
    #     filepath='/Users/kgori/scratch/testcoal',tmpdir='/tmp',gtp_path='/Users/kgori/git/kevin/clustering_project/class_files/',unit_is_pam=True,quiet=True)
