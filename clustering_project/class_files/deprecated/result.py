#!/usr/bin/python
# -*- coding: utf-8 -*-
from collections import defaultdict

class Result(object):

    def update_score(self):
        self.score = sum([rec.tree.score for rec in self.concats])

    def __init__(self, clusters={}):
        self.score = 0
        self.length = 0
        self.concats = []
        self.members = []
        self.names = []
        if clusters:
            self.length = len(clusters)
            for key in sorted(clusters):
                self.concats.append(clusters[key]['concatenation'])
                self.members.append(clusters[key]['members'])
                self.names.append([rec.name for rec in
                                  clusters[key]['members']])
            self.update_score()

    def retrieve_names(self, index):
        return self.names[index]

    def retrieve_concat(self, index):
        return self.concats[index]

    def retrieve_members(self, index):
        return self.members[index]

    def find_mergeable_groups(self, matrix):
        """ 'Matrix' parameter is a numpy array, or python list of lists
            NOT a numpy matrix (indexing works differently)
        """
        nclusters = len(matrix)
        indices = range(nclusters)
        results = defaultdict(list)
        for i in range(nclusters):
            if i in indices:
                for j in range(i+1, nclusters):
                    if matrix[i][j] == 0:
                        results[i].append(j)
                        indices.remove(j)
        results = { k:v for k,v in results.items() if len(v)>0 }
        return (results, matrix)

    def write_ALF_parameters(
        self,
        simulation_name,
        experiment_directory,
        simulation_directory,
        number_of_genes,
        gene_length,
        tree_file,
        outfilename='alf-params.drw',
        ):
        """ Function to write parameter files for ALF """

        import os
        alfsim_parameter_string = \
        '''# This parameter file has been generated by the ALF web service.
# To run the simulation on your own machine, download ALF from http://www.cbrg.ethz.ch/alf
# and call alfsim from the parent directory of this file.

webRequest := true;
uuid := '';

# name of simulation - you may want to change this
mname := {0};

# directories for file storage - you may want to change these
wdir := '{1}/{2}';
dbdir := 'DB/';
dbAncdir := 'DBancestral/';

# time scale for simulation (PAM is default)
unitIsPam := true:

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
'''.format(simulation_name,
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
