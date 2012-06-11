#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import cPickle


def extract(sc, n):
    d = sc.get_clusters()
    key = ('sym', 'ward', n)
    trees = [x.tree.newick for x in d[key].concats]
    lengths = []
    for member_list in d[key].members:
        l = [x.seqlength for x in member_list]
        lengths.append(l)
    total_lengths = [sum(x) for x in lengths]
    out_d = {}
    out_d['trees'] = trees
    out_d['lengths'] = lengths
    out_d['total_lengths'] = total_lengths
    out_d['n'] = n
    return out_d


def write_ALF_parameters(
    simulation_name,
    experiment_directory,
    simulation_directory,
    number_of_genes,
    gene_length,
    tree_file,
    output_filename='alf-params.drw',
    ):
    """ Function to write parameter files for ALF """

    import os
    alfsim_parameter_string = \
        '''# This parameter file has been generated by the ALF web service.
# To run the simulation on your own machine, download ALF from http://www.cbrg.ethz.ch/alf
# and call alfsim from the parent directory of this file.

webRequest := false;
uuid := '';

# name of simulation - you may want to change this
mname := {0};

# directories for file storage - you may want to change these
wdir := '{1}/{2}';
dbdir := 'DB/';
dbAncdir := 'DBancestral/';

# time scale for simulation (PAM is default)
unitIsPam := false:

# parameters concerning the root genome
realseed := false;
protStart := {3};
minGeneLength := {4};
gammaLengthDist := [1,1];
blocksize := 1:

# parameters concerning the substitution models
substModels := [SubstitutionModel('F84', [1, 1], [0.25,0.25,0.25,0.25], false)];
indelModels := [IndelModel(0)];
rateVarModels := [RateVarModel()];
modelAssignments := [1]:
modelSwitchS := [[1]]:
modelSwitchD := [[1]]:

# parameters concerning the species tree
treeType := 'Custom';
treeFile := '{5}';
'''.format(
        simulation_name,
        experiment_directory,
        simulation_directory,
        number_of_genes,
        gene_length,
        tree_file,
        )
    if output_filename:
        if os.path.isfile(output_filename):
            write = \
                raw_input("Output file '{0}' exists, overwrite (y/n)?: ".format(output_filename))
        else:
            write = 'y'
        if write == 'y':
            output = open(output_filename, 'w')
            output.write(alfsim_parameter_string)
            output.close()
    return alfsim_parameter_string


n = int(sys.argv[1])
sc = cPickle.load(file('yeast_jc69.pickle'))
name = '{0}to{1}'.format(n, n + 1)
if not os.path.isdir('./simulations/{0}'.format(name)):
    os.mkdir('./simulations/{0}'.format(name))
simdir = './simulations/{0}'.format(name)

d = extract(sc, n)
lengths = d['lengths']
total_lengths = d['total_lengths']
trees = d['trees']

helper = os.environ['DARWINHELPER']
tmpdir = os.environ['TEMPORARY_DIRECTORY']

parameter_files = []
for i in range(len(total_lengths)):
    treefile = open('{0}/tree{1}.nwk'.format(simdir, i + 1), 'w')
    treefile.write(trees[i])
    treefile.close()
    write_ALF_parameters(
        'alfsim' + name,
        './tmp',
        'alftmp',
        1,
        total_lengths[i],
        '{0}/tree{1}.nwk'.format(simdir, i + 1),
        '{0}/class{1}-params.drw'.format(simdir, i + 1),
        )
    parameter_files.append('{0}/class{1}-params.drw'.format(simdir,
                           i + 1))
    cPickle.dump(lengths[i],file('{0}/lengths{1}.pickle'.format(simdir,i+1),'w'))

