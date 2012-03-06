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
  -m [rf, sym] = matrix: RF or symmetric differences
  -l [single*, complete, average, ward] = linkage method ()
  -cut [float 0 < cut < 1] = cutting point at which to separate clusters (default 0.25)
  -show = plot dendrogram
  -cluster = concatenate sequences based on discovered clusters and write the result to in/clusters
  -norm = normalise matrix before clustering (i.e. divide throughout by maximum value)
''')

#Check arguments
if "-dir" not in args:
    INPUT_DIR = '.'
elif not args['-dir']:
    print "No input directory specified, trying '.' ..."
    INPUT_DIR = ','
else: INPUT_DIR = args['-dir']

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

if "-cluster" not in args:
    make_clusters = False
else: make_clusters = True

if "-norm" not in args:
    normalise = False
else: normalise = True

class RAxML_object(dendropy.Tree):
    lnl = None

taxa = dendropy.TaxonSet()
tree_files = glob.glob( "{0}/trees/besttrees/*".format(INPUT_DIR) )
info_files = glob.glob( "{0}/trees/info/*".format(INPUT_DIR) )
msa_files = glob.glob( "{0}/MSA/*.phy".format(INPUT_DIR) )
likelihoods = [float( re.compile( "(?<=Score of best tree ).+" ).search(open(x).read()).group() ) for x in info_files]
trees = [RAxML_object() for x in tree_files]
num_trees = len(trees)
[trees[i].read_from_path(tree_files[i],'newick',taxon_set=taxa) for i in range(num_trees)]
for i in range(num_trees): 
    trees[i].lnl = likelihoods[i]

np.set_printoptions(precision=2,linewidth=200)

#Set up matrices
matrix = np.zeros( (num_trees,num_trees),dtype='float' )

for i in range(len(trees)):
    for j in range(i+1,len(trees)):
        if matrix_type == 'mix':
            matrix[i][j]=dendropy.treecalc.robinson_foulds_distance(trees[i],trees[j])
            matrix[j][i]=dendropy.treecalc.symmetric_difference(trees[i],trees[j])
        elif matrix_type == 'rf':
            matrix[i][j]=matrix[j][i]=dendropy.treecalc.robinson_foulds_distance(trees[i],trees[j])
        elif matrix_type == 'sym':
            matrix[i][j]=matrix[j][i]=dendropy.treecalc.symmetric_difference(trees[i],trees[j])

if normalise:
    matrix = matrix / np.max(matrix)
    
print matrix
try: 
    Y = squareform(matrix)
    link = linkage(Y, linkage_method)
except: 
    Y = matrix
    link = linkage(Y, linkage_method)

cut = (link[-1][2])*cut
# Plot the clustering via dendrogram function - 'cut' parameter separates the trees into clusters at the given cut-point (set at 1/4 max branch length for convenience)
# T gives a list recording the number of the cluster each gene is assigned to
T = fcluster(link,cut,criterion="distance")
if make_clusters:   
    if not os.path.exists( "{0}/clusters".format(INPUT_DIR)): os.mkdir( "{0}/clusters".format(INPUT_DIR))
    if not os.path.exists( "{0}/clusters/MSA".format(INPUT_DIR)): os.mkdir( "{0}/clusters/MSA".format(INPUT_DIR))
    if not os.path.exists( "{0}/clusters/trees".format(INPUT_DIR)): os.mkdir( "{0}/clusters/trees".format(INPUT_DIR))
    clusters = {}
    for k in range(min(T),max(T)+1):
        clusters[k] = []

    for i in range(len(T)):
        clusters[T[i]].append(get_phylip_file(msa_files[i]))

    for key in clusters.keys():
        seq_list = clusters[key]
        conc = concatenate_alignments(seq_list)
        conc.write_phylip( "{0}/clusters/MSA/cluster{1:0>2}.phy".format(INPUT_DIR,key))

    assignments = zip(tree_files,T)
    assignment_outfile = open( "{0}/clusters/cluster_assignments.txt".format(INPUT_DIR), "w")
    for x in sorted(assignments, key=lambda tup: int(tup[1])): 
        gene = x[0][x[0].rindex("/")+1:x[0].rindex(".")]
        assignment_outfile.write("{0}\t{1}\n".format(gene,x[1]))

if show_plot: 
    dendrogram( link, color_threshold=cut, leaf_font_size=10,leaf_rotation=90,leaf_label_func=lambda leaf: tree_files[leaf][1+tree_files[leaf].rindex('/'):tree_files[leaf].rindex('.')]+"_"+str(T[leaf]),count_sort=True)
    title("{0} linkage of {1} matrix".format(linkage_method,matrix_type))
    axhline(cut,color='grey',ls='dashed')
    xlabel('Gene')
    ylabel('Distance')
    show()
