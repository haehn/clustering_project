#!/usr/bin/env python
import dendropy, glob, numpy, re, sys
from scipy.cluster.hierarchy import single, complete, linkage, dendrogram
from hcluster import squareform
from matplotlib.pyplot import show,title,xlabel,ylabel
from handleArgs import handleArgs
from Bio import Cluster

args = handleArgs(sys.argv,help='''
all_against_all arguments:
  -in = path to data directory (default = '.')
  -m [rf, sym, mix*] = matrix: RF, symmetric differences, or mix of both (upper triangle = RF, lower = symm diff)
  -l [single*, complete, ward] = linkage method ()
''')

#Check arguments
if "-in" not in args:
    INPUT_DIR = '.'
elif not args['-in']:
    print "No input directory specified, trying '.' ..."
    INPUT_DIR = ','
else: INPUT_DIR = args['-in']

if "-m" not in args:
    matrix_type = 'mix'
elif not args['-in']:
    print "No matrix type specified, using 'mix' ..."
    matrix_type = 'mix'
else: matrix_type = args['-m']

if "-l" not in args:
    linkage_method = 'single'
elif not args['-in']:
    print "No linkage method specified, using 'single' ..."
    linkage_method = 'single'
else: linkage_method = args['-l']

class RAxML_object(dendropy.Tree):
    lnl = None

taxa = dendropy.TaxonSet()
tree_files = glob.glob( "{}/trees/besttrees/*".format(INPUT_DIR) )
info_files = glob.glob( "{}/trees/info/*".format(INPUT_DIR) )
likelihoods = [float( re.compile( "(?<=Score of best tree ).+" ).search(open(x).read()).group() ) for x in info_files]
trees = [RAxML_object() for x in tree_files]
num_trees = len(trees)
[trees[i].read_from_path(tree_files[i],'newick',taxon_set=taxa) for i in range(num_trees)]
for i in range(num_trees): 
    trees[i].lnl = likelihoods[i]

numpy.set_printoptions(precision=2,linewidth=200)

#Set up matrices
matrix = numpy.zeros( (num_trees,num_trees),dtype='float' )

for i in range(len(trees)):
    for j in range(i+1,len(trees)):
        if matrix_type == 'mix':
            matrix[i][j]=dendropy.treecalc.robinson_foulds_distance(trees[i],trees[j])
            matrix[j][i]=dendropy.treecalc.symmetric_difference(trees[i],trees[j])
        elif matrix_type == 'rf':
            matrix[i][j]=matrix[j][i]=dendropy.treecalc.robinson_foulds_distance(trees[i],trees[j])
        elif matrix_type == 'sym':
            matrix[i][j]=matrix[j][i]=dendropy.treecalc.symmetric_difference(trees[i],trees[j])

print matrix
try: 
    Y = squareform(matrix)
    link = linkage(Y, linkage_method)
except: 
    link = linkage(matrix, linkage_method)

cut = (link[-1][2])*0.25
dendrogram( link, color_threshold=cut, leaf_font_size=10,leaf_rotation=90,leaf_label_func=lambda leaf: tree_files[leaf][1+tree_files[leaf].rindex('/'):tree_files[leaf].rindex('.')],count_sort=True)
title("{0} linkage of {1} matrix".format(linkage_method,matrix_type))
xlabel('Gene')
ylabel('Distance')
show()
