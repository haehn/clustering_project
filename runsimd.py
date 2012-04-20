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
import rpy2.robjects as rob
import dendropy as dpy
import argparse


# args = handleArgs(sys.argv, help = '''
# runsim.py arguments:
#   -c    =  number of classes
#   -s    =  number of species
#   -nni  =  apply this many NNIs to first class tree to generate subsequent class trees
#   -spr  =  as above, but use SPRs
#   -g    =  number of genes per class, e.g. [10,20,10]
#   -r    =  rate per class (PAM units, default is 200), eg [200,100,400]
#   -dir  =  Output directory to place simulation into (default './')
#   -tmp  =  temporary files from alfsim go here (keep separate when running more than one instance)
# ''')

# #Check arguments: -c should give an integer number of classes, -s an integer number of species, -g a list of number of genes for each class
# # NB This currently doesn't check input type for errors

# C = int(check_args_value( '-c', args, 2  ))
# S = int(check_args_value( '-s', args, 20 ))
# R = int(check_args_value( '-r', args, 250))
# NNI = int(check_args_value( '-nni', args, 0 )) # 0 evaluates to False
# SPR = int(check_args_value('-spr', args, 0 ))
# OUT_DIR = check_args_filepath( '-dir', args, "" )
# TEMP_DIR = check_args_filepath( '-tmp', args, "alftmp" )
# OUT_DIR = os.getcwd()+'/'+OUT_DIR
# TEMP_DIR = os.getcwd()+'/'+TEMP_DIR

# # More complicated cases - no function for this
# if not '-g' in args: 
#     G = [ 10 for each_class in range(C) ]
# elif not args['-g']:
#     print 'Number of genes not specified. Genes set to 10 per class'
#     G = [ 10 for each_class in range(C) ]
# else: 
#     regexp=re.compile('(?<=[\[,])[0-9]+')
#     G = regexp.findall(args['-g']) #Make sure we understand the list given in G
#     if not G: 
#         print "Didn't understand '-g' argument."
#         sys.exit()
#     else: G = [int(x) for x in G]



# #DEBUG ARGS
# #print "Type of G elements: ",type(G[0])
# MSA_path = "{0}/MSA".format(OUT_DIR)
# if not os.path.isdir(MSA_path): os.mkdir(MSA_path)
# tree_path = "{0}/trees".format(OUT_DIR)
# if not os.path.isdir(tree_path): os.mkdir(tree_path)


parser = argparse.ArgumentParser(prog='runsim.py', description='Run ALF simulator to generate different topological classes')
parser.add_argument('-c','--classes', help='Number of classes', type=int, default=2)
parser.add_argument('-s','--species', help='Number of species', type=int, default=20)
parser.add_argument('-g','--genes', help='Number of genes per class', nargs='*', default=[])
parser.add_argument('-r','--rates', help='Tree length per class (PAM units)', nargs='*', default=[])
parser.add_argument('-nni','--nni', help='Number of NNIs', type=int, default=0)
parser.add_argument('-spr','--spr', help='Number of SPRs', type=int, default=0)
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
NNI = args['nni']
SPR = args['spr']
OUT_DIR = args['directory']
TEMP_DIR = args['temp_directory']
MSA_path = fpath("{0}/MSA".format(OUT_DIR))
tree_path = fpath("{0}/trees".format(OUT_DIR))
quiet = args['quiet']
for each in [OUT_DIR, TEMP_DIR, MSA_path, tree_path]:
    if not os.path.isdir(each): os.mkdir(each)

print 'od',OUT_DIR
print 'td',TEMP_DIR
print 'mp',MSA_path
print 'tp',tree_path
print 'nni',NNI
print 'spr',SPR
print G

# Make basetree
write_ALF_parameters(simulation_name='basetree', experiment_directory=TEMP_DIR, \
    simulation_directory='', number_of_genes=1, min_gene_length=250, number_of_species=S, \
    mutation_rate=R[0], indels=True, bd_tree=True, \
    output_filename='{0}/basetree-params.drw'.format(OUT_DIR)) 

run_ALF('{0}/basetree-params.drw'.format(OUT_DIR))
os.rename('{0}/basetree/RealTree.nwk'.format(TEMP_DIR), '{0}/basetree.nwk'.format(tree_path))

#Write parameters
for i in range(C):
    if NNI:
        write_ALF_parameters(simulation_name='class{0}'.format(i+1), experiment_directory=TEMP_DIR, \
            simulation_directory='', number_of_genes=G[i], min_gene_length=250, number_of_species=S, \
            mutation_rate=R[i], indels=True, bd_tree=False, custom_tree='{0}/custom{1}.nwk'.format(tree_path, i+1),\
            output_filename='{0}/class{1}-params.drw'.format(OUT_DIR,i+1))   
    elif SPR:
        write_ALF_parameters(simulation_name='class{0}'.format(i+1), experiment_directory=TEMP_DIR, \
            simulation_directory='', number_of_genes=G[i], min_gene_length=250, number_of_species=S, \
            mutation_rate=R[i], indels=True, bd_tree=False, custom_tree='{0}/custom{1}.nwk'.format(tree_path, i+1),\
            output_filename='{0}/class{1}-params.drw'.format(OUT_DIR,i+1)) 
    else:
        write_ALF_parameters(simulation_name='class{0}'.format(i+1), experiment_directory=TEMP_DIR, \
            simulation_directory='', number_of_genes=G[i], min_gene_length=250, number_of_species=S, \
            mutation_rate=R[i], indels=True, bd_tree=True, output_filename='{0}/class{1}-params.drw'.format(OUT_DIR,i+1))   
    

rob.r('library(phangorn)')
for i in range(C):
    rob.r('t <- read.tree("{0}/basetree.nwk")'.format(tree_path))
    # rob.r('t <- root(t, "SE001",resolve.root=TRUE)')
    rob.r('c <- rNNI(t,{0})'.format(NNI))
    # rob.r('c <- root(c, "SE001",resolve.root=TRUE)')
    rob.r('n <- write.tree(c)')
    n = rob.r['n']
    newick = n[0]
    writer = open("{0}/custom{1}.nwk".format(tree_path,i+1),"w")
    writer.write(newick+';\n')
    writer.close()

for i in range(C):  
    run_ALF( '{0}/class{1}-params.drw'.format(OUT_DIR,i+1) )
    ### DO OUR RELABELLING HERE
    right = dpy.Tree().get_from_stream(open("{0}/custom{1}.nwk".format(tree_path,i+1)),'newick')
    wrong = dpy.Tree().get_from_stream(open("{0}/class{1}/RealTree.nwk".format(TEMP_DIR,i+1)),'newick')
    mismatch_dict = {} # Records 'wrong' label as key and 'right' label as value
    for (n1,n2) in zip(right.preorder_node_iter(),wrong.preorder_node_iter()):
        if n1.is_leaf() and n2.is_leaf():
            mismatch_dict[n2.taxon.label]=n1.taxon.label
            n2.taxon.label=str(n1.taxon)
    wrong.write_to_stream(open('{0}/true{1}.nwk'.format(tree_path,i+1),"w"),"newick")
    #Gather MSA files
    class_path = '{0}/class{1}/MSA'.format(TEMP_DIR,i+1)
    class_msas = glob.glob("{0}/*.fa".format(class_path))
    class_msas = sorted(class_msas, key=lambda name: int(name.split('_')[1]))
    for fi in class_msas: # Give MSA file a better name
        name_elements = fi.split('_')
        gene_number = int(name_elements[1]) + sum(G[:i]) #add on number of genes already named in previous classes
        rename = "{0}/gene{1:0>3}.fas".format(MSA_path, gene_number)
        print fi, rename
        seq_object = get_fasta_file(fi)
        seq_object.headers = [mismatch_dict[x.split('/')[0]] for x in seq_object.headers]
        seq_object.write_fasta(rename)
        os.remove(fi)

    #Gather tree files
    try:
        pam2sps('{0}/true{1}.nwk'.format(tree_path,i+1),"pam2sps",'{0}/true{1}_sps.nwk'.format(tree_path,i+1))
        #os.rename('{0}/class{1}/RealTree.nwk'.format(TEMP_DIR,i+1), '{0}/true{1}.nwk'.format(tree_path,i+1))
    except: pass
                
    #Delete unnecessary simulated filesls
    shutil.rmtree("{0}/class{1}".format(TEMP_DIR,i+1))
shutil.rmtree("{0}/basetree".format(TEMP_DIR,i+1))

true_clustering = open("./{0}/true_clustering".format(MSA_path),"w")
for i in range(len(G)):
    for j in range(G[i]):
        true_clustering.write(str(i+1)+" ")
true_clustering.close()

####### Now we have a labelling problem
# Alf doesn't preserve tip labelling when simulating from a custom user tree, so we need to
# 1: get the correct ordering from the input tree (=true#_sps.nwk)
# 2: apply this labelling to the output tree, and store the mismatches
# 3: relabel the sequence files according to the mismatching

# Step 1
# right = dpy.Tree().get_from_stream("{0}/custom{1}.nwk".format(tree_path,i+1),'newick')
# wrong = dpy.Tree().get_from_stream("{0}/class{1}/RealTree.nwk".format(TEMP_DIR,i+1),'newick')
# mismatch_dict = {} # Records 'wrong' label as key and 'right' label as value
# for (n1,n2) in zip(right.preorder_node_iter(),wrong.preorder_node_iter()):
#     if n1.is_leaf() and n2.is_leaf():
#         mismatch_dict[n2.taxon.label]=n1.taxon.label
#         n2.taxon.label=str(n1.taxon)
# wrong.write_to_stream(open('{0}/true{1}.nwk'.format(tree_path,i+1),"w"),"newick")