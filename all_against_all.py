#!/usr/bin/env python
import dendropy, glob, numpy, sys
from hcluster import single, complete, linkage, dendrogram
from matplotlib.pyplot import show
from handleArgs import handleArgs

args = handleArgs(sys.argv,help='''
all_against_all arguments:
  -in = path to data directory (default = '.')
  -m [rf, sym, mix*] = matrix: RF, symmetric differences, or mix of both (upper triangle = RF, lower = symm diff)
  -l [single*, complete] = linkage method ()
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

trees = glob.glob( "{}/trees/besttrees/*".format(INPUT_DIR) )
#trees = sorted( trees, key = lambda name: int(name.split('/')[-1][4:6]) )
num_trees = len(trees)
taxa = dendropy.TaxonSet()

#Set up matrices
matrix = numpy.zeros( (num_trees,num_trees),dtype='float' )
numpy.set_printoptions(precision=2,linewidth=200)

for i in range(len(trees)):
    tree1 = dendropy.Tree()
    tree1.read_from_path(trees[i],'newick',taxon_set=taxa)
    for j in range(i+1,len(trees)):        
        tree2 = dendropy.Tree()
        tree2.read_from_path(trees[j],'newick',taxon_set=taxa)
        if matrix_type == 'mix':
            matrix[i][j]=dendropy.treecalc.robinson_foulds_distance(tree1,tree2)
            matrix[j][i]=dendropy.treecalc.symmetric_difference(tree1,tree2)
        elif matrix_type == 'rf':
            matrix[i][j]=matrix[j][i]=dendropy.treecalc.robinson_foulds_distance(tree1,tree2)
        elif matrix_type == 'sym':
            matrix[i][j]=matrix[j][i]=dendropy.treecalc.symmetric_difference(tree1,tree2)

print matrix

if linkage_method == 'single': dendrogram( single(matrix),orientation='right' )
elif linkage_method == 'complete': dendrogram( complete(matrix),orientation='right' )

show()















