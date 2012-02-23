#!/usr/bin/env python
import dendropy, glob, numpy
from hcluster import single, complete, linkage, dendrogram
from matplotlib.pyplot import show

trees = glob.glob( "./trees/besttrees/*" )
trees = sorted( trees, key = lambda name: int(name.split('/')[-1][4:6]) )
num_trees = len(trees)
taxa = dendropy.TaxonSet()

mix_matrix = numpy.zeros( (num_trees,num_trees),dtype='float' )
RF_matrix = numpy.zeros( (num_trees,num_trees),dtype='float' )
sym_matrix = numpy.zeros( (num_trees,num_trees),dtype='float' )

for i in range(len(trees)):
    for j in range(i+1,len(trees)):
        tree1 = dendropy.Tree()
        tree2 = dendropy.Tree()
        tree1.read_from_path(trees[i],'newick',taxon_set=taxa)
        tree2.read_from_path(trees[j],'newick',taxon_set=taxa)
        mix_matrix[i][j]=RF_matrix[i][j]=RF_matrix[j][i]=dendropy.treecalc.robinson_foulds_distance(tree1,tree2)
        mix_matrix[j][i]=sym_matrix[i][j]=sym_matrix[j][i]=dendropy.treecalc.symmetric_difference(tree1,tree2)


numpy.set_printoptions(precision=2,linewidth=200)
dendrogram( linkage(mix_matrix),orientation='right' )
#dendrogram( single(mix_matrix),orientation='right' )
#dendrogram( complete(mix_matrix),orientation='right' )
show()


