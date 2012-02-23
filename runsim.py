#!/usr/bin/env python
"""
Quick & dirty script to run alfsim and generate G genes for 20 species, with genes
from C classes (i.e. a mixture of defined histories).
MSAs are put in MSA directory
True trees are put in trees directory
Remaining ALF files are deleted
"""
from ALF_wrapper import *         # functions write_ALF_parameters() and run_ALF()
from sequence_record import *     # classes Unaligned_Sequence_Record and Aligned_Sequence_Record, function get_fasta_file()
from handleArgs import handleArgs
import glob, os, re, shutil, sys

args = handleArgs(sys.argv, help = '''
runsim.py arguments:
  -c  =  number of classes
  -g  =  number of genes (total)
  -s  =  number of species
  -dir = Output directory to place simulation into (default './')
''')

#Check arguments: -c should give an integer number of classes, -s an integer number of species, -g a list of number of genes for each class
# NB This currently doesn't check input type for errors
if not '-c' in args: 
    C = 2
elif not args['-c']: 
    print 'Number of classes not specified. Classes set to 2'
    C = 2
else: C = int(args['-c'])

if not '-s' in args: 
    S = 20
elif not args['-s']:
    print 'Number of species not specified. Species set to 20'
    S = 20
else: S = int(args['-s'])

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

if not '-dir' in args: 
    OUT_DIR = './'
elif not args['-dir']: 
    print 'Output directory not specified. Directory set to current (./)'
    OUT_DIR = './'
else: 
    OUT_DIR = args['-dir']
    if not os.path.isdir(args['-dir']):
        os.mkdir(OUT_DIR)
    else: print "Output Directory already exists, by the way..."

#DEBUG ARGS
print "Type of G elements: ",type(G[0])


MSA_path = "{}/MSA".format(OUT_DIR)
if not os.path.isdir(MSA_path): os.mkdir(MSA_path)
tree_path = "{}/trees".format(OUT_DIR)
if not os.path.isdir(tree_path): os.mkdir(tree_path)

TEMP_DIR = 'scratch'

for i in range(C):
    #Write parameters
    write_ALF_parameters(simulation_name='class{}'.format(i+1), experiment_directory=TEMP_DIR, simulation_directory='', number_of_genes=G[i], min_gene_length=250, number_of_species=S, mutation_rate=250, indels=True, output_filename='{0}/class{1}-params.drw'.format(OUT_DIR,i+1))   

    #Run simulation
    run_ALF('{0}/class{1}-params.drw'.format(OUT_DIR,i+1))

    #Gather MSA files
    class_path = './{0}/class{1}/MSA'.format(TEMP_DIR,i+1)
    class_msas = glob.glob("{0}/*.fa".format(class_path))
    class_msas = sorted(class_msas, key=lambda name: int(name.split('_')[1]))
    for fi in class_msas: # Give MSA file a better name
        name_elements = fi.split('_')
        gene_number = int(name_elements[1]) + sum(G[:i]) #add on number of genes already named in previous classes
        rename = "{0}/MSA_{1:0>2}.fa".format(MSA_path, gene_number)
        print fi, rename
        os.rename(fi, rename)

    #Gather tree files
    os.rename('./{0}/class{1}/RealTree.nwk'.format(TEMP_DIR,i+1), './{0}/true{1}.nwk'.format(tree_path,i+1))
    
    #Delete unnecessary simulated filesls
    shutil.rmtree("./{0}/class{1}".format(TEMP_DIR,i+1))
#END LOOP

#Convert fasta files to phylip files for use with raxml
for fasta in glob.glob("{0}/*.fa".format(MSA_path)):
    phylip = fasta[:fasta.rindex('.')]+'.phy'
    seq_object = get_fasta_file(fasta)
    seq_object.headers = [x.split('/')[0] for x in seq_object.headers]
    seq_object.write_phylip(outfile=phylip, print_to_screen = False)
