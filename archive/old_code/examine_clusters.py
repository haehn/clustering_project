#!/usr/bin/env python
import dendropy, glob, os, re, sys
import numpy as np
from scipy.cluster.hierarchy import single, complete, linkage, dendrogram, fcluster
from hcluster import squareform
from matplotlib.pyplot import show,title,xlabel,ylabel,axhline
from handleArgs import handleArgs
from sequence_record import *

args = handleArgs(sys.argv,help='''
all_against_all arguments:
  -dir = path to data directory (default = '.')
  -clusters = path to cluster trees
  -m [rf, sym] = matrix: RF or symmetric differences
  -l [single*, complete, average, ward] = linkage method ()
  -cut [float 0 < cut < 1] = cutting point at which to separate clusters (default 0.25)
  -show = plot dendrogram
  -norm = normalise matrix before clustering (i.e. divide throughout by maximum value)
''')

#Check arguments
if "-dir" not in args:
    INPUT_DIR = '.'
elif not args['-dir']:
    print "No input directory specified, trying '.' ..."
    INPUT_DIR = ','
else: INPUT_DIR = args['-dir']

if "-clusters" not in args:
    CLUSTER_DIR = '.'
elif not args['-clusters']:
    print "No input directory specified, trying '.' ..."
    CLUSTER_DIR = ','
else: CLUSTER_DIR = args['-clusters']

if "-m" not in args:
    matrix_type = 'rf'
elif not args['-m']:
    print "No matrix type specified, using 'rf' ..."
    matrix_type = 'rf'
else: matrix_type = args['-m']

if "-l" not in args:
    linkage_method = 'single'
elif not args['-l']:
    print "No linkage method specified, using 'single' ..."
    linkage_method = 'single'
else: linkage_method = args['-l']

if "-cut" not in args:
    cut = 0.25
elif not args['-cut']:
    print "Cutting point not specified, using 0.25..."
    cut = 0.25
else: cut = float(args['-cut'])

if "-show" not in args:
    show_plot = False
else: show_plot = True

if "-norm" not in args:
    normalise = False
else: normalise = True

class RAxML_object(dendropy.Tree):
    lnl = None

taxa = dendropy.TaxonSet()
single_tree_files = glob.glob( "{0}/*.nwk".format(INPUT_DIR) )
#single_info_files = glob.glob( "{0}/trees/info/*".format(INPUT_DIR) )
#single_likelihoods = [float( re.compile( "(?<=Score of best tree ).+" ).search(open(x).read()).group() ) for x in single_info_files]
single_trees = [RAxML_object() for x in single_tree_files]
num_genes = len(single_trees)
[single_trees[i].read_from_path(single_tree_files[i],'newick',taxon_set=taxa) for i in range(num_genes)]
#for i in range(num_genes): 
#    single_trees[i].lnl = single_likelihoods[i]
names = [ n[1+n.rindex("/"):n.rindex(".")] for n in single_tree_files]
#print names

cluster_tree_files = glob.glob( "{0}/*.nwk".format(CLUSTER_DIR) )
#cluster_info_files = glob.glob( "{0}/clusters/trees/info/*".format(INPUT_DIR) )
#cluster_likelihoods = [float( re.compile( "(?<=Score of best tree ).+" ).search(open(x).read()).group() ) for x in cluster_info_files]
cluster_trees = [RAxML_object() for x in cluster_tree_files]
num_clusters = len(cluster_trees)
[cluster_trees[i].read_from_path(cluster_tree_files[i],'newick',taxon_set=taxa) for i in range(num_clusters)]
#for i in range(num_clusters): 
#    cluster_trees[i].lnl = cluster_likelihoods[i]

true_tree_files = glob.glob( "{0}../*_sps.nwk".format(INPUT_DIR))
print "{0}../*_sps.nwk".format(INPUT_DIR)
true_trees = [RAxML_object() for x in true_tree_files]
num_true_trees = len(true_trees)
[true_trees[i].read_from_path(true_tree_files[i],'newick',taxon_set=taxa) for i in range(num_true_trees)]


np.set_printoptions(precision=2,linewidth=200)

#Set up matrices
cc_matrix = np.zeros( (num_clusters,num_clusters),dtype='float' ) #Cluster-cluster all against all matrix
comparison_matrix = np.zeros( (num_clusters, num_genes), dtype='float' ) #Cluster-gene comparison matrix
true_concat_matrix = np.zeros( (num_clusters, num_true_trees), dtype='float') #Cluster-true comparison matrix

for i in range(num_clusters):
    for j in range(i+1,num_clusters):
        if matrix_type == 'rf':
            cc_matrix[i][j]=cc_matrix[j][i]=dendropy.treecalc.robinson_foulds_distance(cluster_trees[i],cluster_trees[j])
        elif matrix_type == 'sym':
            cc_matrix[i][j]=cc_matrix[j][i]=dendropy.treecalc.symmetric_difference(cluster_trees[i],cluster_trees[j])
    for k in range(num_genes):
        if matrix_type == 'rf':
            comparison_matrix[i][k]=dendropy.treecalc.robinson_foulds_distance(cluster_trees[i],single_trees[k])
        if matrix_type == 'sym':
            comparison_matrix[i][k]=dendropy.treecalc.symmetric_difference(cluster_trees[i],single_trees[k])
    for l in range(num_true_trees):
        if matrix_type == 'rf':
            true_concat_matrix[i][l]=dendropy.treecalc.robinson_foulds_distance(cluster_trees[i],true_trees[l])
        if matrix_type == 'sym':
            true_concat_matrix[i][l]=dendropy.treecalc.symmetric_difference(cluster_trees[i],true_trees[l])        

if normalise:
    cc_matrix = cc_matrix / np.max(cc_matrix)
    comparison_matrix = comparison_matrix / np.max(comparison_matrix)
    if num_true_trees>0: true_concat_matrix = true_concat_matrix / np.max(true_concat_matrix)
    
#print cc_matrix
#print comparison_matrix



# This dictionary gives us lists of the genes belonging to each cluster (KEY = CLUSTER NUMBER, VALUE = LIST)
clusters = {}

for line in open( "{0}/clusters/cluster_assignments.txt".format(INPUT_DIR) ):
    tup = tuple(line.rstrip().split())
    if len(tup) > 2:
        tup = tup[:2]
    if not tup[1] in clusters:
        clusters[tup[1]] = [ tup[0] ]
    else: clusters[tup[1]].append(tup[0])

#print clusters

zips = []

for i in range(num_clusters):
    z = zip(names,comparison_matrix[i])
    zips.append(z)
print "Average within-groups distance from consensus tree, (min, max):"
for i in range(1,num_clusters+1):
    l = clusters[str(i)]
    print "Cluster {0} ({1} genes)".format(i, len(l))
    for cl in zips:
        scores = []
        for x in cl:
            if x[0] in l:
                scores.append(x[1])
        print "{0:.4f} ({1:.4f}, {2:.4f})".format(np.sum(scores)/len(l), min(scores), max(scores))


if num_true_trees>0:
    print
    print "Correspondence between clusters (rows) and true trees (columns)"
    print true_concat_matrix



try: 
    Y = squareform(cc_matrix)
    link = linkage(Y, linkage_method)
except: 
    Y = cc_matrix
    link = linkage(Y, linkage_method)

cut = (link[-1][2])*cut
# Plot the clustering via dendrogram function - 'cut' parameter separates the trees into clusters at the given cut-point (set at 1/4 max branch length for convenience)
# T gives a list recording the number of the cluster each gene is assigned to
T = fcluster(link,cut,criterion="distance")

if show_plot: 
    dendrogram( link, color_threshold=cut, leaf_font_size=10,leaf_rotation=90,leaf_label_func=lambda leaf: cluster_tree_files[leaf][1+cluster_tree_files[leaf].rindex('/'):cluster_tree_files[leaf].rindex('.')]+"_"+str(T[leaf]),count_sort=True)
    title("{0} linkage of {1} matrix".format(linkage_method,matrix_type))
    axhline(cut,color='grey',ls='dashed')
    xlabel('Gene')
    ylabel('Distance')
    show()
