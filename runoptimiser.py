#!/usr/bin/env python

# Standard libraries
import argparse,glob,os,sys
# External libraries
import numpy as np
import dendropy as dpy
# Own libraries
from utility_functions import *
from inference_functions import *
from sequence_record import *

###################################################################################################

#===========================================#
# Get command line arguments with argparse #
#=========================================#

metrics = ['euc', 'rf', 'sym']
linkages = ['average', 'complete', 'median', 'single', 'ward', 'weighted']
parser = argparse.ArgumentParser(description='Optimise partitioning')
parser.add_argument('-d',   '--directory',          dest='MSA_DIR',         type=fpath, default='./MSA',    metavar='DIRECTORY',    help='Path to MSA files')
parser.add_argument('-t',   '--tree-directory',     dest='TREES_DIR',       type=fpath, default='./tr',     metavar='DIRECTORY',    help='Directory to save output tree files')
parser.add_argument('-c',   '--cluster-directory',  dest='CLUSTER_DIR',     type=fpath, default='./cl',     metavar='DIRECTORY',    help='Directory to save output cluster files')
parser.add_argument('-tmp', '--temp-directory',     dest='TEMP_DIR',        type=fpath, default='.',        metavar='DIRECTORY',    help='Directory to save temp files')
parser.add_argument('-r',   '--results-directory',  dest='RESULTS_DIR',     type=fpath, default='./results',metavar='DIRECTORY',    help='Directory to save output result files')
parser.add_argument('-n',   '--number-of-classes',  dest='nclasses',        type=int,   default=4,                                  help='Number of classes to partition into')
parser.add_argument('-dist','--distance-measure',   dest='distance_metric', type=str,   default='sym',      choices = metrics,      help='Tree distance measure')
parser.add_argument('-l',   '--linkage-type',       dest='linkage_type',    type=str,   default='ward',     choices = linkages,     help='Linkage criterion')
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

for each in [TREES_DIR, CLUSTER_DIR, TEMP_DIR, RESULTS_DIR]:
    if not os.path.isdir(each): os.mkdir(each)

###################################################################################################

#=================================================#
# Read alignments and calculate individual trees #
#===============================================#

fasta_files = get_alignments(MSA_DIR) # we'll need access to the sequence alignments (fasta format)
names = [x[x.rindex("/")+1:x.rindex(".")] for x in fasta_files]

if not quiet: print "OBTAINING GENE TREES"

gene_trees = []

for fasta in fasta_files:
    sanitise_fasta(fasta) # Reorders sequences into alphabetical order, and replaces whitespace with underscores
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
              
    if TREES_DIR:                                       # check for pre-calculated tree result
        try: 
            tree = Inference_Result()
            tree.read_from_file("{0}/{1}.tree".format(TREES_DIR,name))
            if not quiet: print "Found pre-existing tree in output directory: {0}.tree".format( name )
        except IOError:
            if not quiet: print "Running TreeCollection on {0}".format(dv)
            tree = run_treecollection(dv, map_file, labels, guide_tree, name)
            tree.write_to_file( "{0}/{1}.tree".format(TREES_DIR,tree.name) )
    else: 
        if not quiet: print "Running TreeCollection on {0}".format(dv)
        tree = run_treecollection(dv, map_file, labels, guide_tree, name)
    gene_trees.append( tree )

###################################################################################################

#==========================================#
# Cluster trees by topological similarity #
#========================================#

if not quiet: print "CLUSTERING GENE TREES"
matrix = get_distance_matrix(gene_trees,distance_metric,normalise=True)
link = get_linkage(matrix,linkage_type)
if rand:
    clustering = list(np.random.randint(1,int(nclasses)+1,size=len(gene_trees)))
else:
    clustering = cluster_linkage(link,nclasses,criterion='distance')
t = cluster_linkage(link,nclasses,criterion='distance')
if show: showplot(matrix, clustering, link, names, nclasses)
print clustering

if not quiet: print "OBTAINING CLUSTER TREES"
assign_to_clusters_optimiser(fasta_files, clustering, CLUSTER_DIR)
cluster_trees = get_cluster_trees(CLUSTER_DIR)
best_score = sum([float(tr.score) for tr in cluster_trees])

###################################################################################################

#=========================================================#
# Optimisation procedure:
# STARTING POINT
# - clustering is calculated by some means,
#   partitioning genes into n clusters
# - the 'best' tree is calculated for each cluster using
#   info from each gene in the cluster
# - cluster tree scores are summed - objective function
# PROCEDURE
# - a sample of genes are selected
# - for each gene in the sample its individual topology
#   is compared (RF) to the n cluster trees
# - 
#=======================================================#
 
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
resets = 0
while stayed_put < 10:
    if resets == 4:
        print "Reset limit reached ({0})".format(resets)
        break
    if done_worse == 5:
        print "wandered off, resetting..."
        resets+=1
        done_worse = 0
        best_score = global_best[0]
        clustering = global_best[1]
        print i+1, best_score, clustering
    # if stayed_put == 10:
    #     print 'stayed put 10 times, exitting...'
    #     break
    assignments   = assign_to_clusters_optimiser(fasta_files, clustering, CLUSTER_DIR)
    cluster_trees = get_cluster_trees(CLUSTER_DIR, quiet=True)
    best_score    = sum([float(tr.score) for tr in cluster_trees])
    cluster_files = glob.glob("{0}/*".format(CLUSTER_DIR))
    for fi in cluster_files:
        os.remove(fi)

    clustering, best_score=optimise_sample_rf(clustering, gene_trees,\
        cluster_trees, len(gene_trees)/2 ,len(gene_trees)/10, MSA_DIR,\
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
    print i+1, best_score, sc(clustering,t)
    writer.write("{0} {1} {2}\n".format(i+1, best_score, sc(clustering,t)))
    i+=1
writer.close()
#assignments = assign_to_clusters_optimiser(fasta_files, clustering, CLUSTER_DIR)
#cluster_trees = get_cluster_trees(CLUSTER_DIR,quiet=True)
best_score = sum([float(tr.score) for tr in cluster_trees])
print "Best score: {0}".format(global_best[0])
print "Best clustering: {0}".format(global_best[1])

for tr in cluster_trees:
    tr.write_to_file("{0}/{1}.tree".format(CLUSTER_DIR, tr.name))

# pl<-function(df) {
#       plot(df[,1],df[,2],type='l',col='black',lwd=4)
#       par(new=T)
#       plot(df[,1],df[,3],type='l',col='blue',lwd=1,yaxt='n')
#       axis(4)
#       mtext('',side=4,line=3,font=2)
#       par(new=F)
# }
        