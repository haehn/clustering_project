#!/usr/bin/env python
from sequence_record import *
import dendropy, glob, numpy, os, re, sys
from handleArgs import handleArgs

class RAxML_object(dendropy.Tree):
    lnl = None
    def get_lnl(self,raxml_info_file):
        import re
        likelihood = re.compile( "(?<=Score of best tree ).+" ).search(open(raxml_info_file).read()).group()
        self.lnl = float(likelihood)
    def get_raxml_tree(self, raxml_besttree_file, taxon_set=dendropy.TaxonSet()):
        self.read_from_path(raxml_besttree_file, 'newick', taxon_set)
    def call_raxml(self, charset = "aa"):
        import os
        if charset == 'aa': model = "PROTGAMMAWAG"
        elif charset == 'dna': model = "GTRGAMMA"
        pass

args = handleArgs(sys.argv,help='''
all_against_all arguments:
  -in = path to data directory (default = '.')
  -m [rf, sym, mix*] = matrix: RF, symmetric differences, or mix of both (upper triangle = RF, lower = symm diff)
''')

#Check arguments
if "-in" not in args:
    INPUT_DIR = '.'
elif not args['-in']:
    print "No input directory specified, trying '.' ..."
    INPUT_DIR = ','
else: INPUT_DIR = args['-in']

if "-m" not in args:
    matrix_type = 'rf'
elif not args['-in']:
    print "No matrix type specified, using 'rf' ..."
    matrix_type = 'rf'
else: matrix_type = args['-m']

if "-model" not in args:
    model = "PROTGAMMAWAG"
elif not args['-model']:
    print "No sequence type specified, using 'PROTGAMMAWAG' ..."
    model = "PROTGAMMAWAG"
else: model = args['-model']

if "-raxml" not in args:
    run_raxml = False
elif not args['-raxml']:
    print "Using all-against-all raxml"
    run_raxml = True
else: run_raxml = args['-raxml']

#Collect raxml besttree files and info files
#Use these to set up RAxML objects with trees and likelihoods
tree_files = glob.glob( "{0}/trees/besttrees/*".format(INPUT_DIR) )
info_files = glob.glob( "{0}/trees/info/*".format(INPUT_DIR) )
msa_files = glob.glob( "{0}/MSA/*.phy".format(INPUT_DIR))
trees = [RAxML_object() for x in tree_files]
num_trees = len(trees)
taxa = dendropy.TaxonSet()
[trees[i].read_from_path(tree_files[i],'newick',taxon_set=taxa) for i in range(num_trees)]
[trees[i].get_lnl(info_files[i]) for i in range(num_trees)]

numpy.set_printoptions(precision=2,linewidth=200)

#Set up matrices
matrix = numpy.zeros( (num_trees,num_trees),dtype='float' )

max_value = None
min_value = None
for i in range(len(trees)):
    for j in range(i+1,len(trees)):
        if matrix_type == 'rf':
            matrix[i][j]=matrix[j][i]=dendropy.treecalc.robinson_foulds_distance(trees[i],trees[j])
            if not max_value:
                max_value = matrix[i][j]
            else:
                if matrix[i][j] > max_value:
                    max_value = matrix[i][j]
            if not min_value:
                min_value = matrix[i][j]
            else:
                if matrix[i][j] < min_value:
                    min_value = matrix[i][j]

        elif matrix_type == 'sym':
            matrix[i][j]=matrix[j][i]=dendropy.treecalc.symmetric_difference(trees[i],trees[j])
            if not max_value:
                max_value == matrix[i][j]
            else:
                if matrix[i][j] > max_value:
                    max_value = matrix[i][j]
            if not min_value:
                min_value = matrix[i][j]
            else:
                if matrix[i][j] < min_value:
                    min_value = matrix[i][j]

print matrix
print max_value, min_value

if not os.path.exists("{0}/concats".format(INPUT_DIR)): os.mkdir("{0}/concats".format(INPUT_DIR))
if not os.path.exists("{0}/concats/besttrees".format(INPUT_DIR)): os.mkdir("{0}/concats/besttrees".format(INPUT_DIR))
if not os.path.exists("{0}/concats/info".format(INPUT_DIR)): os.mkdir("{0}/concats/info".format(INPUT_DIR))

true1_treeobject = dendropy.Tree()
try: true1_treeobject.read_from_path("{0}/trees/true1_sps.nwk".format(INPUT_DIR),'newick',taxon_set=taxa)
except: true_tree_known = False
true2_treeobject = dendropy.Tree()
try: true2_treeobject.read_from_path("{0}/trees/true2_sps.nwk".format(INPUT_DIR),'newick',taxon_set=taxa)
except: pass

for i in range(len(trees)):
    name_tree_i = tree_files[i][1+tree_files[i].rindex('/'):tree_files[i].rindex('.')]
    for j in range(i+1,len(trees)):
        name_tree_j = tree_files[j][1+tree_files[j].rindex('/'):tree_files[j].rindex('.')]
        name_concat = name_tree_i + "_" + name_tree_j
        concat = concatenate_alignments(get_phylip_file(msa_files[i]),get_phylip_file(msa_files[j]),name=name_concat) # Get these function calls out of the loop!
        concat.write_phylip(outfile="{0}/concats/{1}.phy".format(INPUT_DIR,name_concat))
        sum_lnl = trees[i].lnl + trees[j].lnl
        mutual_distance = matrix[i][j]
        if run_raxml:
            os.system( 'raxml -T 8 -m {2} -s {0}/concats/{1}.phy -n {1} -p 121 && mv RAxML_bestTree.{1} {0}/concats/besttrees/{1}.nwk && mv RAxML_info.{1} {0}/concats/info/{1}.info && rm *.{1} '.format(INPUT_DIR,name_concat, model) )
        	#print 'raxml -T 8 -m {2} -s {0}/concats/{1}.phy -n {1} -p 121 && mv RAxML_bestTree.{1} {0}/concats/besttrees/{1}.nwk && mv RAxML_info.{1} {0}/concats/info/{1}.info && rm *.{1} '.format(INPUT_DIR,name_concat, model)
        tree_i_object = trees[i]
        tree_j_object = trees[j]
        concat_tree_object = dendropy.Tree()
        concat_tree_object.read_from_path("{0}/concats/besttrees/{1}.nwk".format(INPUT_DIR,name_concat),'newick',taxon_set=taxa)
        concat_lnl = float(re.compile( "(?<=Score of best tree ).+" ).search(open("{0}/concats/info/{1}.info".format(INPUT_DIR,name_concat)).read()).group())
        dist_tree1_concat = dendropy.treecalc.robinson_foulds_distance(trees[i],concat_tree_object)
        dist_tree2_concat = dendropy.treecalc.robinson_foulds_distance(trees[j],concat_tree_object)
        if true_tree_known:
            dist_true1_tree_i = dendropy.treecalc.robinson_foulds_distance(true1_treeobject,trees[i])
            dist_true2_tree_i = dendropy.treecalc.robinson_foulds_distance(true2_treeobject,trees[i])
            dist_true1_tree_j = dendropy.treecalc.robinson_foulds_distance(true1_treeobject,trees[j])
            dist_true2_tree_j = dendropy.treecalc.robinson_foulds_distance(true2_treeobject,trees[j])
            dist_true1_concat = dendropy.treecalc.robinson_foulds_distance(true1_treeobject,concat_tree_object)
            dist_true2_concat = dendropy.treecalc.robinson_foulds_distance(true2_treeobject,concat_tree_object)
            print "{0}{1:15.4f}{2:15.4f}{3:15.4f}{4:10.4f}{5:10.4f}{6:10.4f}{7:10.4f}{8:10.4f}{9:10.4f}{10:10.4f}".format( name_concat, sum_lnl, concat_lnl, concat_lnl - sum_lnl, mutual_distance, dist_true1_tree_i, dist_true2_tree_i, dist_true1_tree_j, dist_true2_tree_j, dist_true1_concat, dist_true2_concat )
        else: 
            print "{0}{1:15.4f}{2:15.4f}{3:15.4f}{4:10.4f}{5:10.4f}{6:10.4f}".format( name_concat, sum_lnl, concat_lnl, concat_lnl - sum_lnl, mutual_distance, dist_tree1_concat, dist_tree2_concat )
        