#!/usr/bin/env python

import glob,os,sys
import numpy as np
import dendropy as dpy
from utility_functions import *
from inference_functions import *
from sequence_record import *

##################################################################################################

command_line_args = handleArgs(sys.argv,help='''
runml.py arguments:
         -program      = raxml, phyml or treecollection
         -dna          = invokes GTRGAMMA model in RAxML (default is PROTGAMMAWAG)
         -dir          = input directory
         -v            = verbose (include spawned process stdout) 
         -out          = directory to write gene trees to
         -nclasses     = number of classes to cluster into
         --cluster-out = path to write clusters to
         --use-results = path to results directory if gene trees have already been calculated
         --use-cluster = path to cluster trees if already calculated
''')

#Check command line arguments

if not '-program' in command_line_args:
    print "No program chosen"
    sys.exit(1)
elif not command_line_args['-program']:
    print "No program chosen"
    sys.exit(1)
else:
    program = command_line_args["-program"]
    # else:
    #     print "Unrecognised program chosen - use 'raxml', 'phyml' or 'treecollection'"
    #     sys.exit(1)

if not '-dir' in command_line_args:
    INPUT_DIR = '.'
elif not command_line_args['-dir']:
    print "No input directory specified, trying '.'"
    INPUT_DIR = '.'
else: 
    INPUT_DIR = command_line_args['-dir']
    while INPUT_DIR.endswith("/"): INPUT_DIR = INPUT_DIR[:-1] # Trim off trailing "/"s

if '--use-results' in command_line_args:
    try: 
        GENE_TREES_DIR = command_line_args["--use-results"]
        while GENE_TREES_DIR.endswith("/"): GENE_TREES_DIR = GENE_TREES_DIR[:-1] # Trim off trailing "/"s
        results = True
    except KeyError:
        results = False
else: results = False

if '-out' in command_line_args:
    try: 
        OUT_DIR = command_line_args["-out"]
        while OUT_DIR.endswith("/"): OUT_DIR = OUT_DIR[:-1] # Trim off trailing "/"s
        write_genetrees = True
    except KeyError:
        write_genetrees = False
else: write_genetrees = False

if '-nclasses' in command_line_args:
    try: nclasses = command_line_args['-nclasses']
    except KeyError: nclasses = 4
else: nclasses = 4 

if '--write-clusters' in command_line_args:
    try: 
        CLUSTER_DIR = command_line_args["--write-clusters"]
        while CLUSTER_DIR.endswith("/"): CLUSTER_DIR = CLUSTER_DIR[:-1] # Trim off trailing "/"s
    except KeyError:
        CLUSTER_DIR = None
else: CLUSTER_DIR = None

if '-dna' in command_line_args: 
    dna = True
else:
    dna = False 

if '-v' in command_line_args:
    verbose = True
else: 
    verbose = False

##################################################################################################

fasta_files = get_alignments(INPUT_DIR) # we'll need access to the sequence alignments (fasta format)
phylip_files = [x[:x.rindex(".")]+".phy" for x in fasta_files]

if results: # avoid wasting time by checking for inference results in trees directory
    gene_trees = list(get_gene_trees(GENE_TREES_DIR))

else: # if no results already, make trees
    gene_trees = []

    for fasta in fasta_files:
        prefix = fasta[:fasta.rindex(".")]
        phylip = prefix+".phy"
        dv = prefix+"_dv.txt"
        if dna:
            datatype = 'DNA'
        else: datatype = 'AA'  
        if not os.path.isfile(phylip):
                print "Making phylip file for ", fasta
                populate_phylip_from_fasta(fasta) 
        if not os.path.isfile(dv):
                print "making TreeCollection input files for ", fasta
                populate_dv_from_fasta(fasta,datatype)

        if program == 'raxml':            
            if dna: model = "GTRGAMMA"
            else: model = "PROTGAMMAWAG"
            alignment = phylip
            name = prefix[prefix.rindex("/")+1:]
            print "Running RAxML on {0}".format(phylip)
            gene_trees.append( run_raxml( model, alignment, name, nthreads=8 ) )
        
        elif program == 'phyml':
            if dna:
                model = 'GTR'
                datatype = 'nt'
            else: 
                model = 'WAG'
                datatype = 'aa' 
            alignment = phylip
            name = prefix[prefix.rindex("/")+1:]
            print "Running PhyML on {0}".format(phylip)
            gene_trees.append( run_phyml( model, alignment, name, datatype) )
            
        elif program == 'treecollection':          
            map_file = prefix+"_map.txt"
            labels = prefix+"_labels.txt"
            guide_tree = prefix+"_guidetree.nwk"
            name = prefix[prefix.rindex("/")+1:]
            print "Running TreeCollection on {0}".format(dv)
            gene_trees.append( run_treecollection(dv, map_file, labels, guide_tree, name) )
   
if write_genetrees:
    if not os.path.isdir(OUT_DIR): os.mkdir(OUT_DIR)
    for tree in gene_trees:
        tree.write_to_file("{0}/{1}.trobj".format(OUT_DIR,tree.name))

""" Fourth: clustering """
matrix = get_distance_matrix(gene_trees)
link = get_linkage(matrix,"single")
clustering = cluster_linkage(link,nclasses)
print clustering
assignments = assign_to_clusters(phylip_files, clustering, CLUSTER_DIR)

""" Fifth: generate some stats """
cluster_files = get_alignments(CLUSTER_DIR)
cluster_trees = []

for cluster in cluster_files:
    prefix = cluster[:cluster.rindex(".")]
    phylip = prefix+".phy"
    dv = prefix+"_dv.txt"

    if program == 'raxml':
        if dna: model = "GTRGAMMA"
        else: model = "PROTGAMMAWAG"
        alignment = phylip
        name = prefix[prefix.rindex("/")+1:]
        print "Running RAxML on {0}".format(phylip)
        cluster_trees.append( run_raxml( model, alignment, name, nthreads=8 ) )
    
    elif program == 'phyml':
        if dna:
            model = 'GTR'
            datatype = 'nt'
        else: 
            model = 'WAG'
            datatype = 'aa' 
        alignment = phylip
        name = prefix[prefix.rindex("/")+1:]
        print "Running PhyML on {0}".format(phylip)
        cluster_trees.append( run_phyml( model, alignment, name, datatype) )
        
    elif program == 'treecollection':
        if dna:
            datatype = 'DNA'
        else: datatype = 'AA'
        dv = prefix+"_dv.txt"
        map_file = prefix+"_map.txt"
        labels = prefix+"_labels.txt"
        guide_tree = prefix+"_guidetree.nwk"
        name = prefix[prefix.rindex("/")+1:]
        print "Running TreeCollection on {0}".format(dv)
        cluster_trees.append( run_treecollection(dv, map_file, labels, guide_tree, name) )

if write_genetrees:
    if not os.path.isdir(OUT_DIR): os.mkdir(OUT_DIR)
    for tree in cluster_trees:
        tree.write_to_file("{0}/{1}.trobj".format(OUT_DIR,tree.name))

#  - MESSY CODE

gtdict = {}
for tree in gene_trees:
    gtdict[tree.name]=tree

for key in assignments:
    print "Cluster",key
    trees = []
    for gene in assignments[key]:
        ir_object = gtdict[gene] # Inference_Result object (i.e. a tree object)
        trees.append(ir_object)
    
    # Get average within cluster difference
    dm = get_distance_matrix(trees)
    triu = np.triu_indices(len(dm),1)
    print "Mean of all_against_all RF (topology) comparison: ",np.mean(dm[triu])

    # Get concat tree distance from true tree
    tr_tree = False
    taxa = dpy.TaxonSet()
    cl_tree = dpy.Tree()
    cl_tree.read_from_string(cluster_trees[key-1].tree,'newick',taxon_set=taxa)
    true_trees = glob.glob( "{0}/../trees/true*_sps.nwk".format(INPUT_DIR))
    trcomps = []
    for tr_tree in true_trees:
        tr_tree_dpy = dpy.Tree()
        try: tr_tree_dpy.read_from_string(open(tr_tree).read(),'newick',taxon_set=taxa)
        except IOError: 
            tr_tree = False
            break
        trcomps.append(tr_tree_dpy.symmetric_difference(cl_tree))
       
    if tr_tree: print "Cluster vs True RF (topology) comparison", min(trcomps), trcomps

    # Get mean distance from concat tree
    dpytrees = [dpy.Tree() for tree in trees]
    for i in range(len(trees)):
        dpytrees[i].read_from_string(trees[i].tree,'newick',taxon_set=taxa)
    meancomps = []
    for dpytree in dpytrees:
        meancomps.append(dpytree.symmetric_difference(cl_tree))
    print "Average RF (topology) distance from concat tree: ",np.mean(meancomps)
