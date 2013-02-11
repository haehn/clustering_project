#!/usr/bin/env python

from textwrap import dedent


class Params(object):

    """ Class for writing parameter files for use with alfsim """

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
        """ "One-word" models don't need any extra parameters. Codon: CPAM, ECM,
        ECMu AA:    GCB,  JTT, LG,  WAG """

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
        gain_rate=0.00001,
        gain_model='ZIPF',
        gain_params=[1.821],
        max_gain_length=10,
        loss_rate=0.00001,
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
        """ Models rate variation among sites Only uses gamma distribution (ALF
        allows more options) """

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
        outfile_name = '{0}/{1}.drw'.format(self.outfile_path, self.name)
        writer = open(outfile_name, 'w')
        writer.write(str(self))
        writer.flush()
        writer.close()

        return outfile_name
