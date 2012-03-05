#!/usr/bin/env python
"""
Quick write to disk of cluster dendrograms of all-against-all tree distance matrix for exploratory analysis
"""

import dendropy, glob, numpy, sys
from scipy.cluster.hierarchy import single, complete, linkage, dendrogram, fcluster
from hcluster import squareform
from matplotlib.pyplot import clf,savefig,show,title,xlabel,ylabel,axhline
from handleArgs import handleArgs
from Bio import Cluster
from math import sqrt, trunc

args = handleArgs(sys.argv,help='''
all_against_all arguments:
  -in = path to data directory (default = '.')
  -out = save prefix
  -cut [float 0 < cut < 1] = cutting point at which to separate clusters (default 0.25)
''')

#Check arguments
if "-in" not in args:
    INPUT_DIR = '.'
elif not args['-in']:
    print "No input directory specified, trying '.' ..."
    INPUT_DIR = ','
else: INPUT_DIR = args['-in']

if "-out" not in args:
    save_prefix = 'out'
    args['-out'] = 'out'
elif not args['-out']:
    print "No save prefix specified, using 'out' ..."
    save_prefix = 'out'
else: save_prefix = args['-out']

if "-n" not in args:
    nclusters = None
elif not args['-out']:
    print "No save prefix specified, using 'out' ..."
    nclusters= None
else: nclusters = int(args['-n'])

if "-cut" not in args:
    cut_proportion = 0.25
elif not args['-in']:
    print "Cutting point not specified, using 0.25..."
    cut_proportion = 0.25
else: cut_proportion = float(args['-cut'])

taxa = dendropy.TaxonSet()
tree_files = glob.glob( "{}/trees/besttrees/*".format(INPUT_DIR) )
trees = [dendropy.Tree() for x in tree_files]
num_trees = len(trees)
[trees[i].read_from_path(tree_files[i],'newick',taxon_set=taxa) for i in range(num_trees)]
if not nclusters: nclusters = int(sqrt(num_trees*0.5))
if num_trees > 80: font_size = 4
else: font_size = 8

numpy.set_printoptions(precision=2,linewidth=200)

#Set up matrices
rf_matrix = numpy.zeros( (num_trees,num_trees),dtype='float' )
sym_matrix = numpy.zeros( (num_trees,num_trees),dtype='float' )

for i in range(len(trees)):
    for j in range(i+1,len(trees)):          
        rf_matrix[i][j]=rf_matrix[j][i]=dendropy.treecalc.robinson_foulds_distance(trees[i],trees[j])
        sym_matrix[i][j]=sym_matrix[j][i]=dendropy.treecalc.symmetric_difference(trees[i],trees[j])  

linkages = ['single','complete','average','weighted','ward']
matrices = [rf_matrix, sym_matrix]
matrix_names = ['rf','sym']

for x in range(len(linkages)):
    for y in range(len(matrices)):
        filename = "{0}{1}_{2}_{3}.pdf".format(INPUT_DIR,save_prefix,linkages[x],matrix_names[y])
        #k = Cluster.kmedoids(matrices[y],npass=10,nclusters=nclusters)
        try: 
            Y = squareform(matrices[y])
            link = linkage(Y, linkages[x])
        except:
            Y = matrices[y]
            link = linkage(Y, linkages[x])
        cut = (link[-1][2])*cut_proportion
        T = fcluster(link,cut,criterion="distance")
        dendrogram( link, color_threshold=cut, leaf_font_size=font_size, leaf_rotation=90,leaf_label_func=lambda leaf: tree_files[leaf][1+tree_files[leaf].rindex('/'):tree_files[leaf].rindex('.')]+"_"+str(T[leaf]),count_sort=True)
        title("{0} linkage of {1} matrix".format(linkages[x],matrix_names[y]))
        axhline(cut,color='grey',ls='dashed')
        xlabel('Gene')
        ylabel('Distance')
        savefig(filename,format='pdf',dpi=1600)
        clf()       
