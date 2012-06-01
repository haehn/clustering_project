#!/usr/bin/env python

import argparse,glob,os,sys
import numpy as np
import dendropy as dpy
from utility_functions import *
from inference_functions import *
from sequence_record import *

###############################################################################

#============================================#
# Get command line arguments with argparse  #
#==========================================#

metrics  = ['euc', 'rf', 'sym', 'geo']
linkages = ['average', 'complete', 'median', 'single', 'ward', 'weighted']
programs = ["raxml", "phyml", "treecollection"]
parser = argparse.ArgumentParser(description='Optimise partitioning')
parser.add_argument('-d',   '--directory',          dest='MSA_DIR',         type=fpath, default='./MSA',    metavar='DIRECTORY',    help='Path to MSA files')
parser.add_argument('-t',   '--tree-directory',     dest='TREES_DIR',       type=fpath, default='./tr',     metavar='DIRECTORY',    help='Directory to save output tree files')
parser.add_argument('-c',   '--cluster-directory',  dest='CLUSTER_DIR',     type=fpath, default='./cl',     metavar='DIRECTORY',    help='Directory to save output cluster files')
parser.add_argument('-tmp', '--temp-directory',     dest='TEMP_DIR',        type=fpath, default='.',        metavar='DIRECTORY',    help='Directory to save temp files')
parser.add_argument('-r',   '--results-directory',  dest='RESULTS_DIR',     type=fpath, default='./results',metavar='DIRECTORY',    help='Directory to save output result files')
parser.add_argument('-n',   '--number-of-classes',  dest='nclasses',        type=int,   default=4,                                  help='Number of classes to partition into')
parser.add_argument('-dist','--distance-measure',   dest='distance_metric', type=str,   default='sym',      choices = metrics,      help='Tree distance measure')
parser.add_argument('-l',   '--linkage-type',       dest='linkage_type',    type=str,   default='ward',     choices = linkages,     help='Linkage criterion')
parser.add_argument('-program',                     dest='program',         type=str,   default='ward',     choices = programs,     help='Program to use')
parser.add_argument('-dna',                         dest='dna',             action='store_true',                                    help='Nucleotide data')
parser.add_argument('-v',   '--verbose',            dest='verbose',         action='store_true',                                    help='Lots of printing to stdout')
parser.add_argument('-q',   '--quiet',              dest='quiet',           action='store_true',                                    help='Not much printing to stdout')
parser.add_argument('-s',   '--show',               dest='show',            action='store_true',                                    help='Plot dendrograms with pylab')
parser.add_argument('-rand','--randomise',          dest='rand',            action='store_true',                                    help='Start with random partitioning')
args = vars(parser.parse_args())

MSA_DIR         = args['MSA_DIR']
TREES_DIR       = args['TREES_DIR']
CLUSTER_DIR     = args['CLUSTER_DIR']
TEMP_DIR        = args['TEMP_DIR']
RESULTS_DIR     = args['RESULTS_DIR']
nclasses        = args['nclasses']
distance_metric = args['distance_metric']
linkage_type    = args['linkage_type']
dna             = args['dna']
verbose         = args['verbose']
quiet           = args['quiet']
show            = args['show']
rand            = args['rand']
program         = args['program']

for each in [TREES_DIR, CLUSTER_DIR, TEMP_DIR, RESULTS_DIR]:
    if not os.path.isdir(each): os.mkdir(each)

#=#

###############################################################################

fasta_files = get_alignments(MSA_DIR) # we'll need access to the sequence alignments (fasta format)
phylip_files = [x[:x.rindex(".")]+".phy" for x in fasta_files]
names = [x[x.rindex("/")+1:x.rindex(".")] for x in fasta_files]

if not quiet: print "OBTAINING GENE TREES"

gene_trees = []
for fasta in fasta_files:
    sanitise_fasta(fasta)
    prefix = fasta[:fasta.rindex(".")]
    phylip = prefix+".phy"
    dv = prefix+"_dv.txt"
    if dna:
        datatype = 'DNA'
    else: datatype = 'AA'  
    if not os.path.isfile(phylip):
        if not quiet: print "Writing {0} in phylip format...".format(fasta)
        populate_phylip_from_fasta(fasta) 
    if not os.path.isfile(dv):
        if not quiet: print "Making TreeCollection input files for {0}...".format(fasta)
        populate_dv_from_fasta(fasta,datatype)

    if program == 'raxml':                                      # Refactor this code into function
        if dna: model = "GTRGAMMA"
        else: model = "PROTGAMMAWAG"
        alignment = phylip
        name = prefix[prefix.rindex("/")+1:]
        if TREES_DIR:                                       # check for pre-calculated tree result
            try: 
                tree = Inference_Result()
                tree.read_from_file("{0}/{1}.tree".format(TREES_DIR,name))
                if not quiet: print "Found pre-existing tree in output directory"
                assert tree.program.lower() == program                    # 
            except IOError:
                if not quiet: print "Running RAxML on {0}".format(phylip)
                tree = run_raxml( model, alignment, name, nthreads=8 )
                tree.write_to_file( "{0}/{1}.tree".format(TREES_DIR,tree.name) )
            except AssertionError:
                if not quiet: print "This tree was generated by", tree.program
        else: 
            if not quiet: print "Running RAxML on {0}".format(phylip)
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
        if TREES_DIR:                                       # check for pre-calculated tree result
            try: 
                tree = Inference_Result()
                tree.read_from_file("{0}/{1}.tree".format(TREES_DIR,name))
                if not quiet: print "Found pre-existing tree in output directory"
                assert tree.program.lower() == program                    # 
            except IOError:
                if not quiet: print "Running PhyML on {0}".format(phylip)
                tree = run_phyml( model, alignment, name, datatype)
                tree.write_to_file( "{0}/{1}.tree".format(TREES_DIR,tree.name) )
            except AssertionError:
                if not quiet: print "This tree was generated by", tree.program
        else: 
            if not quiet: print "Running PhyML on {0}".format(phylip)
            tree = run_phyml( model, alignment, name, datatype)
        gene_trees.append( tree )

        
    elif program == 'treecollection':          
        map_file = prefix+"_map.txt"
        labels = prefix+"_labels.txt"
        guide_tree = prefix+"_guidetree.nwk"
        name = prefix[prefix.rindex("/")+1:]
        if TREES_DIR:                                       # check for pre-calculated tree result
            try: 
                tree = Inference_Result()
                tree.read_from_file("{0}/{1}.tree".format(TREES_DIR,name))
                if not quiet: print "Found pre-existing tree in output directory: {0}.tree".format( name )
                assert tree.program.lower() == program                    # 
            except IOError:
                if not quiet: print "Running TreeCollection on {0}".format(dv)
                tree = run_treecollection(dv, map_file, labels, guide_tree, name)
                tree.write_to_file( "{0}/{1}.tree".format(TREES_DIR,tree.name) )
            except AssertionError:
                if not quiet: print "This tree was generated by", tree.program
        else: 
            if not quiet: print "Running TreeCollection on {0}".format(dv)
            tree = run_treecollection(dv, map_file, labels, guide_tree, name)
        gene_trees.append( tree )


""" Fourth: clustering """
if not quiet: print "CLUSTERING GENE TREES"

matrix = get_distance_matrix(gene_trees,distance_metric,normalise=True, tmpdir=TEMP_DIR)
link = get_linkage(matrix,linkage_type)
if rand:
    clustering = order(list(np.random.randint(1,int(nclasses)+1,size=len(gene_trees))))
else:
    clustering = order(cluster_linkage(link,nclasses,criterion='distance'))
if show: showplot(matrix, clustering, link, names, nclasses)
assignments = assign_to_clusters(phylip_files, clustering, CLUSTER_DIR)

if not quiet: print "OBTAINING CLUSTER TREES"
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
        if TREES_DIR:                                       # check for pre-calculated tree result
            try: 
                tree = Inference_Result()
                tree.read_from_file("{0}/{1}.tree".format(TREES_DIR,name))
                if not quiet: print "Found pre-existing tree in output directory"
                assert tree.program.lower() == program                    # 
            except IOError:
                if not quiet: print "Running RAxML on {0}".format(phylip)
                tree = run_raxml( model, alignment, name, nthreads=8 )
                tree.write_to_file( "{0}/{1}.tree".format(TREES_DIR,tree.name) )
            except AssertionError:
                if not quiet: print "This tree was generated by", tree.program
        else: 
            if not quiet: print "Running RAxML on {0}".format(phylip)
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
        if TREES_DIR:                                       # check for pre-calculated tree result
            try: 
                tree = Inference_Result()
                tree.read_from_file("{0}/{1}.tree".format(TREES_DIR,name))
                if not quiet: print "Found pre-existing tree in output directory"
                assert tree.program.lower() == program                    # 
            except IOError:
                if not quiet: print "Running PhyML on {0}".format(phylip)
                tree = run_phyml( model, alignment, name, datatype)
                tree.write_to_file( "{0}/{1}.tree".format(TREES_DIR,tree.name) )
            except AssertionError:
                if not quiet: print "This tree was generated by", tree.program
        else: 
            if not quiet: print "Running PhyML on {0}".format(phylip)
            tree = run_phyml( model, alignment, name, datatype)
        cluster_trees.append( tree )

        
    elif program == 'treecollection':          
        map_file = prefix+"_map.txt"
        labels = prefix+"_labels.txt"
        guide_tree = prefix+"_guidetree.nwk"
        name = prefix[prefix.rindex("/")+1:]
        if TREES_DIR:                                       # check for pre-calculated tree result
            try: 
                tree = Inference_Result()
                tree.read_from_file("{0}/{1}.tree".format(TREES_DIR,name))
                if not quiet: print "Found pre-existing tree in output directory"
                assert tree.program.lower() == program                    # 
            except IOError:
                if not quiet: print "Running TreeCollection on {0}".format(dv)
                tree = run_treecollection(dv, map_file, labels, guide_tree, name)
                tree.write_to_file( "{0}/{1}.tree".format(TREES_DIR,tree.name) )
            except AssertionError:
                if not quiet: print "This tree was generated by", tree.program
        else: 
            if not quiet: print "Running TreeCollection on {0}".format(dv)
            tree = run_treecollection(dv, map_file, labels, guide_tree, name)
        cluster_trees.append( tree )

""" Check for any matching cluster trees """
cl_matrix = get_distance_matrix(cluster_trees,"sym")
if not quiet: print cl_matrix
groups = calc_distinct_groups(cl_matrix)

""" Fifth: generate some stats """ #  - MESSY CODE
print "STATS"
if os.path.isfile("{0}/true_clustering".format(MSA_DIR)): 
    print "Linkage: {0}\nMetric: {1}".format(linkage_type,distance_metric)
    #if max(clustering) < 8: print "Misclassifications:", sc([str(x) for x in list(clustering)],open("{0}/true_clustering".format(MSA_DIR)).read().split())
    VI = variation_of_information(clustering,open("{0}/true_clustering".format(MSA_DIR)).read().split())
    print "Variation of Information: ", VI
    # print order([str(x) for x in list(clustering)]),open("{0}/true_clustering".format(MSA_DIR)).read().split()
print "Method: ", program
gtdict = {}                  # Use tree name to pull out tree object
for tree in gene_trees:
    gtdict[tree.name]=tree

index = 0
for key in assignments:
    print "Cluster {0}\nSize: {1}\n{2}".format(key,len(assignments[key]), assignments[key])
    trees = []
    for gene in assignments[key]:
        ir_object = gtdict[gene] # Inference_Result object (i.e. a tree object)
        trees.append(ir_object)
    
    # Get average within cluster difference
    dm = get_distance_matrix(trees)
    if len(dm)>1: 
        triu = np.triu_indices(len(dm),1)
        print "Mean, variance of all against all RF (topology) comparison: ",np.mean(dm[triu]), np.var(dm[triu])

    # Get concat tree distance from true tree
    tr_tree = False
    taxa = dpy.TaxonSet()
    cl_tree = dpy.Tree()
    cl_tree.read_from_string(cluster_trees[key-1].tree,'newick',taxon_set=taxa)
    true_trees = glob.glob( "{0}/../trees/true*_sps.nwk".format(MSA_DIR))
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
    print "Mean, variance of RF (topology) distance from concat tree: {0} {1}".format(np.mean(meancomps),np.var(meancomps))
    print meancomps
    
    # meancomps_tr = []
    # for tree in trees:
    #     meancomps_tr.append(dpy.Tree().get_from_string(tree.tree,'newick').symmetric_difference(dpy.Tree().get_from_stream(open(true_trees[index]),'newick')))
    # print "Mean, variance of RF (topology) distance from true tree: {0} {1}\n".format(np.mean(meancomps_tr),np.var(meancomps_tr))
    # index += 1
    # print meancomps_tr

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
 

writer = open("{0}/result_{1}_{2}.txt".format(RESULTS_DIR, distance_metric, linkage_type),'a')
writer.write("{0} {1}\n".format(VI, final_score))
writer.flush()
writer.close()

        