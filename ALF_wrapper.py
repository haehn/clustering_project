#!/usr/bin/env python

def write_ALF_parameters(simulation_name, experiment_directory, 
    simulation_directory,  number_of_genes=20, min_gene_length=100, 
    number_of_species=20, mutation_rate=200, indels=True, 
    output_filename="alf-params.drw"):
    """ Function to write parameter files for ALF """
    import os
    alfsim_parameter_string = '''# This parameter file has been generated by the ALF web service.
# To run the simulation on your own machine, download ALF from http://www.cbrg.ethz.ch/alf
# and call alfsim from the parent directory of this file.

webRequest := false;
uuid := '';

# name of simulation - you may want to change this
mname := {0};

# directories for file storage - you may want to change these
wdir := './{1}/{2}';
dbdir := 'DB/';
dbAncdir := 'DBancestral/';

# time scale for simulation (PAM is default)
unitIsPam := true:

# parameters concerning the root genome
realseed := false;
protStart := {3};
minGeneLength := {4};
gammaLengthDist := [1, 1];
blocksize := 1:

# parameters concerning the species tree
treeType := 'BDTree';
birthRate := 0.01;
deathRate := 0.001;
NSpecies := {5};
ultrametric := false;
mutRate := {6};
scaleTree := false;

# parameters concerning the substitution models
substModels := [SubstitutionModel('WAG')];
indelModels := [IndelModel(0.0001, ZIPF, [1.821], 10)];
rateVarModels := [RateVarModel()];
modelAssignments := [1]:
modelSwitchS := [[1]]:
modelSwitchD := [[1]]: '''.format(simulation_name, experiment_directory, simulation_directory, number_of_genes, min_gene_length, number_of_species, mutation_rate)
    if not indels:
        alfsim_parameter_string = alfsim_parameter_string.replace("indelModels := [IndelModel(0.0001, ZIPF, [1.821], 10)];\n","")
    if output_filename:
        if os.path.isfile(output_filename):
            write=raw_input("Output file '{0}' exists, overwrite (y/n)?: ".format(output_filename))
        else: write = 'y'
        if write == 'y':
            output = open(output_filename, 'w')
            output.write(alfsim_parameter_string)
            output.close()
    return alfsim_parameter_string


def run_ALF(parameters):
    """ Function to run ALF from parameter file """
    import os
    if os.path.isfile(parameters):
        os.system("alfsim {0}".format(parameters))
    else: print "Can't file file '{0}'".format(parameters)
