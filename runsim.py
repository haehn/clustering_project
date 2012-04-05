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
import glob, os, re, shutil, sys

args = handleArgs(sys.argv, help = '''
runsim.py arguments:
  -c   =  number of classes
  -s   =  number of species
  -g   =  number of genes per class, e.g. [10,20,10]
  -r   =  rate per class (PAM units, default is 200), eg [200,100,400]
  -dir =  Output directory to place simulation into (default './')
  -tmp =  temporary files from alfsim go here (keep separate when running more than one instance)
''')

#Check arguments: -c should give an integer number of classes, -s an integer number of species, -g a list of number of genes for each class
# NB This currently doesn't check input type for errors

C = int(check_args_value( '-c', args, 2  ))
S = int(check_args_value( '-s', args, 20 ))
OUT_DIR = check_args_filepath( '-dir', args, "./" )
TEMP_DIR = check_args_filepath( '-tmp', args, "alftmp" )
print TEMP_DIR

# More complicated cases - no function for this
if not '-g' in args: 
    G = [ 10 for each_class in range(C) ]
elif not args['-g']:
    print 'Number of genes not specified. Genes set to 10 per class'
    G = [ 10 for each_class in range(C) ]
else: 
    regexp=re.compile('(?<=[\[,])[0-9]+')
    G = regexp.findall(args['-g']) #Make sure we understand the list given in G
    if not G: 
        print "Didn't understand '-g' argument."
        sys.exit()
    else: G = [int(x) for x in G]

if not '-r' in args: 
    R = [ 200 for each_class in range(C) ]
elif not args['-r']:
    print 'Mutation rate not specified. Rate set to 200 for each class'
    R = [ 200 for each_class in range(C) ]
else: 
    R = regexp.findall(args['-r']) #Make sure we understand the list given in R
    if not R: 
        print "Didn't understand '-r' argument."
        sys.exit()
    else: R = [int(x) for x in R]

#DEBUG ARGS
#print "Type of G elements: ",type(G[0])

MSA_path = "{0}/MSA".format(OUT_DIR)
if not os.path.isdir(MSA_path): os.mkdir(MSA_path)
tree_path = "{0}/trees".format(OUT_DIR)
if not os.path.isdir(tree_path): os.mkdir(tree_path)

for i in range(C):
    #Write parameters
    write_ALF_parameters(simulation_name='class{0}'.format(i+1), experiment_directory=TEMP_DIR, simulation_directory='', number_of_genes=G[i], min_gene_length=250, number_of_species=S, mutation_rate=R[i], indels=True, output_filename='{0}/class{1}-params.drw'.format(OUT_DIR,i+1))   

    #Run simulation
    run_ALF('{0}/class{1}-params.drw'.format(OUT_DIR,i+1))

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
        
    #Delete unnecessary simulated filesls
    #shutil.rmtree("./{0}/class{1}".format(TEMP_DIR,i+1))
    shutil.rmtree("./{0}/".format(TEMP_DIR))
#END LOOP