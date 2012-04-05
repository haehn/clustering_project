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
         -dna       = invokes GTRGAMMA model in RAxML (default is PROTGAMMAWAG)
         -dir       = input directory
         -v         = verbose (include spawned process stdout)
         -q         = quiet mode: suppress (most) print statements
         -trees     = directory to write gene trees to
         -clusters  = path to write clusters to
         -nclasses  = number of classes to cluster into
         -distance  = sym, rf or euc - tree distance measures
         -rand      = generate random cluster assignment
''')

#Check command line arguments

INPUT_DIR               = check_args_filepath( "-dir", command_line_args, "./MSA" )
OUT_DIR                 = check_args_filepath( "-trees", command_line_args, None )
CLUSTER_DIR             = check_args_filepath( "-clusters", command_line_args, None ) 
TMP_DIR                 = check_args_filepath( "-temp", command_line_args, "./" )
RESULTS_DIR             = check_args_filepath( "-results", command_line_args, None )
dna                     = check_args_bool('-dna', command_line_args)
verbose                 = check_args_bool('-v', command_line_args)
show                    = check_args_bool('-show', command_line_args)
randomise_clustering    = check_args_bool('-rand', command_line_args) 
quiet                   = check_args_bool('-q', command_line_args)
nclasses                = int(check_args_value('-nclasses', command_line_args, 4 ))
matrix_type             = check_args_value('-m', command_line_args, 'sym')
linkage_type            = check_args_value('-l', command_line_args, 'ward')
if quiet: show = False
try:    
    assert nclasses >= 1
    assert matrix_type in ["rf","euc","sym"]
    assert linkage_type in ["average","single","complete","ward","weighted","centroid","median"]
except: sys.exit(1)

##################################################################################################

fasta_files = get_alignments(INPUT_DIR) # we'll need access to the sequence alignments (fasta format)
names = [x[x.rindex("/")+1:x.rindex(".")] for x in fasta_files]

if not quiet: print "OBTAINING GENE TREES"

gene_trees = []

for fasta in fasta_files:
    sanitise_fasta(fasta)
    prefix      = fasta[:fasta.rindex(".")]
    name        = prefix[prefix.rindex("/")+1:]
    dv          = prefix+"_dv.txt"
    map_file    = prefix+"_map.txt"
    labels      = prefix+"_labels.txt"
    guide_tree  = prefix+"_guidetree.nwk"
    if dna:
        datatype    = 'DNA'
    else: 
        datatype    = 'AA'

    if not os.path.isfile(dv):
        if not quiet: print "Making TreeCollection input files for {0}...".format(fasta)
        populate_dv_from_fasta(fasta,datatype)       
              
    if OUT_DIR:                                       # check for pre-calculated tree result
        try: 
            tree = Inference_Result()
            tree.read_from_file("{0}/{1}.tree".format(OUT_DIR,name))
            if not quiet: print "Found pre-existing tree in output directory: {0}.tree".format( name )
        except IOError:
            if not quiet: print "Running TreeCollection on {0}".format(dv)
            tree = run_treecollection(dv, map_file, labels, guide_tree, name)
            tree.write_to_file( "{0}/{1}.tree".format(OUT_DIR,tree.name) )
    else: 
        if not quiet: print "Running TreeCollection on {0}".format(dv)
        tree = run_treecollection(dv, map_file, labels, guide_tree, name)
    gene_trees.append( tree )


""" Fourth: clustering """
if not quiet: print "CLUSTERING GENE TREES"
matrix = get_distance_matrix(gene_trees,matrix_type,normalise=True)
link = get_linkage(matrix,linkage_type)
if randomise_clustering:
    clustering = list(np.random.randint(1,int(nclasses)+1,size=len(gene_trees)))
else:
    clustering = cluster_linkage(link,nclasses,criterion='distance')
if show: showplot(matrix, clustering, link, names, nclasses)
print clustering



if not quiet: print "OBTAINING CLUSTER TREES"
assign_to_clusters_optimiser(fasta_files, clustering, CLUSTER_DIR)
cluster_trees = get_cluster_trees(CLUSTER_DIR)
best_score = sum([float(tr.score) for tr in cluster_trees])



"""
By now there is a clustering list assigning each gene to a cluster, and a score associated with this clustering
Can try to improve this clustering with optimise_clustering() function.

"""

 
suffix = 0    
while os.path.isfile("{0}/result{1}.txt".format(RESULTS_DIR, suffix)):
    suffix += 1
writer = open("{0}/result{1}.txt".format(RESULTS_DIR, suffix),'w')

print "OPTIMISING"
print '0', best_score, clustering


global_best = (best_score, clustering)
done_worse = 0
stayed_put = 0
i=0
while stayed_put < 10:
    if done_worse == 5:
        print "wandered off, resetting..."
        done_worse = 0
        best_score = global_best[0]
        clustering = global_best[1]
        print i+1, best_score, clustering
    # if stayed_put == 10:
    #     print 'stayed put 10 times, exitting...'
    #     break
    assignments   = assign_to_clusters_optimiser(fasta_files, clustering, CLUSTER_DIR)
    print assignments
    sys.exit()
    cluster_trees = get_cluster_trees(CLUSTER_DIR, quiet=True)
    best_score    = sum([float(tr.score) for tr in cluster_trees])
    cluster_files = glob.glob("{0}/*".format(CLUSTER_DIR))
    for fi in cluster_files:
        os.remove(fi)

    clustering, best_score=optimise_sample_rf(clustering, gene_trees,\
        cluster_trees, len(gene_trees)/2, len(gene_trees)/2, INPUT_DIR,\
        CLUSTER_DIR, fasta_files, best_score, greedy=False)

    if best_score < global_best[0]:
        global_best = (best_score, clustering)
        stayed_put = 0
        done_worse = 0
    elif best_score == global_best[0]:
        stayed_put += 1
        done_worse = 0
    else:
        stayed_put = 0
        done_worse += 1
    print i+1, best_score, clustering
    writer.write("{0} {1} {2}\n".format(i+1, best_score, clustering))
    i+=1
writer.close()
#assignments = assign_to_clusters_optimiser(fasta_files, clustering, CLUSTER_DIR)
#cluster_trees = get_cluster_trees(CLUSTER_DIR,quiet=True)
best_score = sum([float(tr.score) for tr in cluster_trees])
print "Best score: {0}".format(global_best[0])
print "Best clustering: {0}".format(global_best[1])

for tr in cluster_trees:
    tr.write_to_file("{0}/{1}.tree".format(CLUSTER_DIR, tr.name))
        