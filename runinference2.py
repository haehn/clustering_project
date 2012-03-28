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
         -program   = raxml, phyml or treecollection
         -dna       = invokes GTRGAMMA model in RAxML (default is PROTGAMMAWAG)
         -dir       = input directory
         -v         = verbose (include spawned process stdout) 
         -trees     = directory to write gene trees to
         -clusters  = path to write clusters to
         -nclasses  = number of classes to cluster into
         -distance  = sym, rf or euc - tree distance measures
         -rand      = generate random cluster assignment
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

if '-trees' in command_line_args:
    try: 
        OUT_DIR = command_line_args["-trees"]
        while OUT_DIR.endswith("/"): OUT_DIR = OUT_DIR[:-1] # Trim off trailing "/"s
        if not os.path.isdir(OUT_DIR): os.mkdir(OUT_DIR)
        write_genetrees = True
    except KeyError:
        write_genetrees = False
else: write_genetrees = False

if '-nclasses' in command_line_args:
    try: 
        nclasses = int(command_line_args['-nclasses'])
        assert nclasses > 0
    except KeyError: nclasses = 4
    except AssertionError: nclasses = 1
    except ValueError: 
        print "{0} is not a valid choice for number of classes".format(command_line_args['-nclasses'])
        sys.exit(0)
else: nclasses = 4 

if '-clusters' in command_line_args:
    try: 
        CLUSTER_DIR = command_line_args["-clusters"]
        while CLUSTER_DIR.endswith("/"): CLUSTER_DIR = CLUSTER_DIR[:-1] # Trim off trailing "/"s
    except KeyError:
        CLUSTER_DIR = None
else: CLUSTER_DIR = None

if "-m" in command_line_args:
    try:
        matrix_type = command_line_args["-m"]
        if not matrix_type in ["rf","euc","sym"]:
            matrix_type = "sym"
    except KeyError:
        matrix_type = "sym"
else: matrix_type = "sym"

if "-l" in command_line_args:
    try:
        linkage_type = command_line_args["-l"]
        if not linkage_type in ["average","single","complete", "ward", "weighted", "centroid", "median"]:
            linkage_type = "ward"
    except KeyError:
        linkage_type = "ward"
else: linkage_type = "ward"

if '-dna' in command_line_args: 
    dna = True
else:
    dna = False 

if '-v' in command_line_args:
    verbose = True
else: 
    verbose = False

if "-show" in command_line_args:
    show = True
else: 
    show = False

if "-rand" in command_line_args:
    randomise_clustering = True
else: 
    randomise_clustering = False

##################################################################################################

fasta_files = get_alignments(INPUT_DIR) # we'll need access to the sequence alignments (fasta format)
phylip_files = [x[:x.rindex(".")]+".phy" for x in fasta_files]
names = [x[x.rindex("/")+1:x.rindex(".")] for x in fasta_files]

print "OBTAINING GENE TREES"

gene_trees = []

for fasta in fasta_files:
    prefix = fasta[:fasta.rindex(".")]
    phylip = prefix+".phy"
    dv = prefix+"_dv.txt"
    if dna:
        datatype = 'DNA'
    else: datatype = 'AA'  
    if not os.path.isfile(phylip):
            print "Writing {0} in phylip format...".format(fasta)
            populate_phylip_from_fasta(fasta) 
    if not os.path.isfile(dv):
            print "Making TreeCollection input files for {0}...".format(fasta)
            populate_dv_from_fasta(fasta,datatype)

    if program == 'raxml':                                      # Refactor this code into function
        if dna: model = "GTRGAMMA"
        else: model = "PROTGAMMAWAG"
        alignment = phylip
        name = prefix[prefix.rindex("/")+1:]
        if write_genetrees:                                       # check for pre-calculated tree result
            try: 
                tree = Inference_Result()
                tree.read_from_file("{0}/{1}.tree".format(OUT_DIR,name))
                print "Found pre-existing tree in output directory"
                assert tree.program.lower() == program                    # 
            except IOError:
                print "Running RAxML on {0}".format(phylip)
                tree = run_raxml( model, alignment, name, nthreads=8 )
                tree.write_to_file( "{0}/{1}.tree".format(OUT_DIR,tree.name) )
            except AssertionError:
                print "This tree was generated by", tree.program
        else: 
            print "Running RAxML on {0}".format(phylip)
            tree = run_raxml( model, alignment, name, nthreads=8 )
        gene_trees.append( tree )
        
    
    elif program == 'phyml':
        if dna:
            model = 'GTR'
            datatype = 'nt'
        else: 
            model = 'WAG'
            datatype = 'aa' 
        alignment = phylip
        name = prefix[prefix.rindex("/")+1:]
        if write_genetrees:                                       # check for pre-calculated tree result
            try: 
                tree = Inference_Result()
                tree.read_from_file("{0}/{1}.tree".format(OUT_DIR,name))
                print "Found pre-existing tree in output directory"
                assert tree.program.lower() == program                    # 
            except IOError:
                print "Running PhyML on {0}".format(phylip)
                tree = run_phyml( model, alignment, name, datatype)
                tree.write_to_file( "{0}/{1}.tree".format(OUT_DIR,tree.name) )
            except AssertionError:
                print "This tree was generated by", tree.program
        else: 
            print "Running PhyML on {0}".format(phylip)
            tree = run_phyml( model, alignment, name, datatype)
        gene_trees.append( tree )

        
    elif program == 'treecollection':          
        map_file = prefix+"_map.txt"
        labels = prefix+"_labels.txt"
        guide_tree = prefix+"_guidetree.nwk"
        name = prefix[prefix.rindex("/")+1:]
        if write_genetrees:                                       # check for pre-calculated tree result
            try: 
                tree = Inference_Result()
                tree.read_from_file("{0}/{1}.tree".format(OUT_DIR,name))
                print "Found pre-existing tree in output directory: {0}.tree".format( name )
                assert tree.program.lower() == program                    # 
            except IOError:
                print "Running TreeCollection on {0}".format(dv)
                tree = run_treecollection(dv, map_file, labels, guide_tree, name)
                tree.write_to_file( "{0}/{1}.tree".format(OUT_DIR,tree.name) )
            except AssertionError:
                print "This tree was generated by", tree.program
        else: 
            print "Running TreeCollection on {0}".format(dv)
            tree = run_treecollection(dv, map_file, labels, guide_tree, name)
        gene_trees.append( tree )


""" Fourth: clustering """
print "CLUSTERING GENE TREES"

matrix = get_distance_matrix(gene_trees,matrix_type,normalise=True)
link = get_linkage(matrix,linkage_type)
if randomise_clustering:
    clustering = list(np.random.randint(1,int(nclasses)+1,size=len(gene_trees)))
else:
    clustering = cluster_linkage(link,nclasses,criterion='distance')
if show: showplot(matrix, clustering, link, names, nclasses)
assignments = assign_to_clusters(phylip_files, clustering, CLUSTER_DIR)

print "OBTAINING CLUSTER TREES"
cluster_files = get_alignments(CLUSTER_DIR)
cluster_names = [ x[x.rindex("/")+1:x.rindex(".")] for x in cluster_files]
cluster_trees = []

for cluster in cluster_files:
    prefix = cluster[:cluster.rindex(".")]
    phylip = prefix+".phy"
    dv = prefix+"_dv.txt"

    if program == 'raxml':                                      # Refactor this code into function
        if dna: model = "GTRGAMMA"
        else: model = "PROTGAMMAWAG"
        alignment = phylip
        name = prefix[prefix.rindex("/")+1:]
        if write_genetrees:                                       # check for pre-calculated tree result
            try: 
                tree = Inference_Result()
                tree.read_from_file("{0}/{1}.tree".format(OUT_DIR,name))
                print "Found pre-existing tree in output directory"
                assert tree.program.lower() == program                    # 
            except IOError:
                print "Running RAxML on {0}".format(phylip)
                tree = run_raxml( model, alignment, name, nthreads=8 )
                tree.write_to_file( "{0}/{1}.tree".format(OUT_DIR,tree.name) )
            except AssertionError:
                print "This tree was generated by", tree.program
        else: 
            print "Running RAxML on {0}".format(phylip)
            tree = run_raxml( model, alignment, name, nthreads=8 )
        cluster_trees.append( tree )
        
    
    elif program == 'phyml':
        if dna:
            model = 'GTR'
            datatype = 'nt'
        else: 
            model = 'WAG'
            datatype = 'aa' 
        alignment = phylip
        name = prefix[prefix.rindex("/")+1:]
        if write_genetrees:                                       # check for pre-calculated tree result
            try: 
                tree = Inference_Result()
                tree.read_from_file("{0}/{1}.tree".format(OUT_DIR,name))
                print "Found pre-existing tree in output directory"
                assert tree.program.lower() == program                    # 
            except IOError:
                print "Running PhyML on {0}".format(phylip)
                tree = run_phyml( model, alignment, name, datatype)
                tree.write_to_file( "{0}/{1}.tree".format(OUT_DIR,tree.name) )
            except AssertionError:
                print "This tree was generated by", tree.program
        else: 
            print "Running PhyML on {0}".format(phylip)
            tree = run_phyml( model, alignment, name, datatype)
        cluster_trees.append( tree )

        
    elif program == 'treecollection':          
        map_file = prefix+"_map.txt"
        labels = prefix+"_labels.txt"
        guide_tree = prefix+"_guidetree.nwk"
        name = prefix[prefix.rindex("/")+1:]
        if write_genetrees:                                       # check for pre-calculated tree result
            try: 
                tree = Inference_Result()
                tree.read_from_file("{0}/{1}.tree".format(OUT_DIR,name))
                print "Found pre-existing tree in output directory"
                assert tree.program.lower() == program                    # 
            except IOError:
                print "Running TreeCollection on {0}".format(dv)
                tree = run_treecollection(dv, map_file, labels, guide_tree, name)
                tree.write_to_file( "{0}/{1}.tree".format(OUT_DIR,tree.name) )
            except AssertionError:
                print "This tree was generated by", tree.program
        else: 
            print "Running TreeCollection on {0}".format(dv)
            tree = run_treecollection(dv, map_file, labels, guide_tree, name)
        cluster_trees.append( tree )

""" Check for any matching cluster trees """
cl_matrix = get_distance_matrix(cluster_trees,"sym",normalise=True)
print cl_matrix
groups = calc_distinct_groups(cl_matrix)

""" Fifth: generate some stats """ #  - MESSY CODE
print "STATS"
gtdict = {}
for tree in gene_trees:
    gtdict[tree.name]=tree

for key in assignments:
    print "Cluster {0}\nSize: {1}\n{2}".format(key,len(assignments[key]), assignments[key])
    trees = []
    for gene in assignments[key]:
        ir_object = gtdict[gene] # Inference_Result object (i.e. a tree object)
        trees.append(ir_object)
    
    # Get average within cluster difference
    dm = get_distance_matrix(trees)
    triu = np.triu_indices(len(dm),1)
    print "Mean of all against all RF (topology) comparison: ",np.mean(dm[triu])

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
       
    if tr_tree: print "Cluster vs True RF (topology) comparison:", min(trcomps), trcomps

    # Get mean distance from concat tree
    dpytrees = [dpy.Tree() for tree in trees]
    for i in range(len(trees)):
        dpytrees[i].read_from_string(trees[i].tree,'newick',taxon_set=taxa)
    meancomps = []
    for dpytree in dpytrees:
        meancomps.append(dpytree.symmetric_difference(cl_tree))
    print "Average RF (topology) distance from concat tree: {0}\n".format(np.mean(meancomps))

# Get sum of likelihood scores
individual_score = sum([float(tr.score) for tr in gene_trees])
final_score = sum([float(tr.score) for tr in cluster_trees])
print "Of {0} clusters, {1} are distinct".format(len (cl_matrix), groups)
print "Individual score: {0}".format(individual_score)
print "Final score: {0}".format(final_score)
if len(cl_matrix) > 1:    
    cluster_link = get_linkage(cl_matrix,"ward")
    cluster_clustering = cluster_linkage(cluster_link, groups, criterion="distance")
    if show: showplot(cl_matrix,cluster_clustering,cluster_link,cluster_names, groups)
