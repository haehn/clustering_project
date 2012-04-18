#!/usr/bin/env python
"""
Script to run alfsim and generate G genes for 20 species, with genes from C classes (i.e. a mixture of defined histories).
MSAs are put in MSA directory
True trees are put in trees directory
Remaining ALF files are deleted
"""
from ALF_wrapper import *         # functions write_ALF_parameters() and run_ALF()
from sequence_record import *     # classes Unaligned_Sequence_Record and Aligned_Sequence_Record, function get_fasta_file()
from utility_functions import * 
import argparse, glob, os, re, shutil, sys

####################################################################################################

#=============================#
# get command line arguments #
#===========================#

parser = argparse.ArgumentParser(prog='runsim.py', description='Run ALF simulator to generate different topological classes')
parser.add_argument('-c','--classes', help='Number of classes', type=int, default=2)
parser.add_argument('-s','--species', help='Number of species', type=int, default=20)
parser.add_argument('-g','--genes', help='Number of genes per class', nargs='*', default=[])
parser.add_argument('-r','--rates', help='Tree length per class (PAM units)', nargs='*', default=[])
parser.add_argument('-d','-dir','--directory', help='Base output directory', type=fpath, default='.')
parser.add_argument('-tmp','--temp-directory', help='Directory to use for temp files', type=fpath, default='./alftmp')
parser.add_argument('-q','--quiet',dest='quiet',action='store_true',help='Not much printing to stdout')
args = vars(parser.parse_args())

# Do some tidying up
while len(args['genes']) < args['classes']: # Make sure each class has a number of genes associated with it
    args['genes'].append(10)
while len(args['rates']) < args['classes']: # Make sure each class has a rate
    args['rates'].append(200)

# Shorter variable names (CAPS = command line args)
C = args['classes']
S = args['species']
G = [int(x) for x in args['genes']]
R = [int(y) for y in args['rates']]
OUT_DIR = args['directory']
TEMP_DIR = args['temp_directory']
MSA_path = fpath("{0}/MSA".format(OUT_DIR))
tree_path = fpath("{0}/trees".format(OUT_DIR))
quiet = args['quiet']
for each in [OUT_DIR, TEMP_DIR, MSA_path, tree_path]:
    if not os.path.isdir(each): os.mkdir(each)

####################################################################################################

#============#
# MAIN LOOP #
#==========#

for i in range(C):
    #Write parameters
    write_ALF_parameters(simulation_name='class{0}'.format(i+1), experiment_directory=TEMP_DIR, simulation_directory='', number_of_genes=G[i], min_gene_length=250, number_of_species=S, mutation_rate=R[i], indels=True, output_filename='{0}/class{1}-params.drw'.format(OUT_DIR,i+1))   

    #Run simulation
    run_ALF('{0}/class{1}-params.drw'.format(OUT_DIR,i+1),quiet)

    #Gather MSA files
    class_path = './{0}/class{1}/MSA'.format(TEMP_DIR,i+1)
    class_msas = glob.glob("{0}/*.fa".format(class_path))
    class_msas = sorted(class_msas, key=lambda name: int(name.split('_')[1]))
    for fi in class_msas: # Give MSA file a better name
        name_elements = fi.split('_')
        gene_number = int(name_elements[1]) + sum(G[:i]) #add on number of genes already named in previous classes
        rename = "{0}/gene{1:0>3}.fas".format(MSA_path, gene_number)
        print fi, rename
        seq_object = get_fasta_file(fi)
        seq_object.headers = [x.split('/')[0] for x in seq_object.headers]
        seq_object.write_fasta(rename)
        os.remove(fi)

    #Gather tree files
    pam2sps('./{0}/class{1}/RealTree.nwk'.format(TEMP_DIR,i+1),"pam2sps",'./{0}/true{1}_sps.nwk'.format(tree_path,i+1))
    os.rename('./{0}/class{1}/RealTree.nwk'.format(TEMP_DIR,i+1), './{0}/true{1}.nwk'.format(tree_path,i+1))
        
    #Delete temporary files
    shutil.rmtree("./{0}/".format(TEMP_DIR))

####################################################################################################

true_clustering = open("./{0}/true_clustering".format(MSA_path),"w")
for i in range(len(G)):
    for j in range(G[i]):
        true_clustering.write(str(i+1)+" ")
true_clustering.close()

