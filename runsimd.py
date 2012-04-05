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


args = handleArgs(sys.argv, help = '''
runsim.py arguments:
  -c    =  number of classes
  -s    =  number of species
  -nni  =  apply this many NNIs to first class tree to generate subsequent class trees
  -spr  =  as above, but use SPRs
  -g    =  number of genes per class, e.g. [10,20,10]
  -r    =  rate per class (PAM units, default is 200), eg [200,100,400]
  -dir  =  Output directory to place simulation into (default './')
  -tmp  =  temporary files from alfsim go here (keep separate when running more than one instance)
''')

#Check arguments: -c should give an integer number of classes, -s an integer number of species, -g a list of number of genes for each class
# NB This currently doesn't check input type for errors

C = int(check_args_value( '-c', args, 2  ))
S = int(check_args_value( '-s', args, 20 ))
R = int(check_args_value( '-r', args, 250))
NNI = int(check_args_value( '-nni', args, 0 )) # 0 evaluates to False
SPR = int(check_args_value('-spr', args, 0 ))
OUT_DIR = check_args_filepath( '-dir', args, "" )
TEMP_DIR = check_args_filepath( '-tmp', args, "alftmp" )
OUT_DIR = os.getcwd()+'/'+OUT_DIR
TEMP_DIR = os.getcwd()+'/'+TEMP_DIR

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



#DEBUG ARGS
#print "Type of G elements: ",type(G[0])
MSA_path = "{0}/MSA".format(OUT_DIR)
if not os.path.isdir(MSA_path): os.mkdir(MSA_path)
tree_path = "{0}/trees".format(OUT_DIR)
if not os.path.isdir(tree_path): os.mkdir(tree_path)


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
    mutation_rate=R, indels=True, bd_tree=True, \
    output_filename='{0}/basetree-params.drw'.format(OUT_DIR)) 

run_ALF('{0}/basetree-params.drw'.format(OUT_DIR))
os.rename('{0}/basetree/RealTree.nwk'.format(TEMP_DIR), '{0}/basetree.nwk'.format(tree_path))

#Write parameters
for i in range(C):
    if NNI:
        write_ALF_parameters(simulation_name='class{0}'.format(i+1), experiment_directory=TEMP_DIR, \
            simulation_directory='', number_of_genes=G[i], min_gene_length=250, number_of_species=S, \
            mutation_rate=R, indels=True, bd_tree=False, custom_tree='{0}/true{1}_sps.nwk'.format(tree_path, i+1),\
            output_filename='{0}/class{1}-params.drw'.format(OUT_DIR,i+1))   
    elif SPR:
        write_ALF_parameters(simulation_name='class{0}'.format(i+1), experiment_directory=TEMP_DIR, \
            simulation_directory='', number_of_genes=G[i], min_gene_length=250, number_of_species=S, \
            mutation_rate=R, indels=True, bd_tree=False, custom_tree='{0}/true{1}_sps.nwk'.format(tree_path, i+1),\
            output_filename='{0}/class{1}-params.drw'.format(OUT_DIR,i+1)) 
    else:
        write_ALF_parameters(simulation_name='class{0}'.format(i+1), experiment_directory=TEMP_DIR, \
            simulation_directory='', number_of_genes=G[i], min_gene_length=250, number_of_species=S, \
            mutation_rate=R, indels=True, bd_tree=True, output_filename='{0}/class{1}-params.drw'.format(OUT_DIR,i+1))   
    
# Do NNI/SPR rearrangements on basetree using R phangorn library
# rob.r('library(phangorn)')
# for i in range(C):
#     cmd = 'n <- write.tree(rNNI(read.tree("{0}/basetree.nwk"),{1}))'.format(tree_path,NNI)
#     rob.r(cmd)
#     n = rob.r['n']
#     tr = dpy.Tree().get_from_string(n[0],'newick')
#     for node in tr.preorder_node_iter():
#         if str(node.taxon) == 'SE001':
#             tr.reroot_at_node(node.parent_node)
#     newick = tr.as_newick_string()
#     writer = open("{0}/true{1}_sps.nwk".format(tree_path,i+1),"w")
#     writer.write(newick+';\n')
#     writer.close()

rob.r('library(phangorn)')
for i in range(C):
    print "AAAAA"
    print "BBBBBB"
    print "CCCCCCC"
    rob.r('t <- read.tree("{0}/basetree.nwk")'.format(tree_path))
    rob.r('t <- root(t, "SE001",resolve.root=TRUE)')
    rob.r('c <- rNNI(t,{0})'.format(NNI))
    rob.r('c <- root(c, "SE001",resolve.root=TRUE)')
    rob.r('n <- write.tree(c)')
    n = rob.r['n']
    newick = n[0]
    writer = open("{0}/true{1}_sps.nwk".format(tree_path,i+1),"w")
    writer.write(newick+';\n')
    writer.close()

for i in range(C):  
    run_ALF( '{0}/class{1}-params.drw'.format(OUT_DIR,i+1) )
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
        seq_object.headers = [x.split('/')[0] for x in seq_object.headers]
        seq_object.write_fasta(rename)
        os.remove(fi)

    #Gather tree files
    try:
        pam2sps('{0}/class{1}/RealTree.nwk'.format(TEMP_DIR,i+1),"pam2sps",'{0}/true{1}_sps.nwk'.format(tree_path,i+1))
        os.rename('{0}/class{1}/RealTree.nwk'.format(TEMP_DIR,i+1), '{0}/true{1}.nwk'.format(tree_path,i+1))
    except: pass
                
    #Delete unnecessary simulated filesls
    shutil.rmtree("{0}/class{1}".format(TEMP_DIR,i+1))
    shutil.rmtree("{0}/basetree".format(TEMP_DIR,i+1))
