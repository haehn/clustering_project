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
  -dir = path to tree files (default = '.')
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

if "-inv" not in args: #invert symmetric difference matrix s.t. each distance is subtracted from the maximum difference,
    invert = False     #which is 2n - 6 where n = number of leaves on the tree
else: invert = True

# DEFINITIONS

class RAxML_object(dendropy.Tree):
    lnl = None

def get_distance_matrix(trees, matrix_type="sym", invert=False, normalise=False):
    import numpy as np
    """ all pairwise distances between trees in list: 'trees' """
    num_trees = len(trees)
    matrix = np.zeros( (num_trees,num_trees),dtype='float' )
    for i in range(num_trees):
        for j in range(i+1,num_trees):
            if matrix_type == 'rf':
                matrix[i][j]=matrix[j][i]=trees[i].robinson_foulds_distance(trees[j])
            elif matrix_type == 'sym':
                matrix[i][j]=matrix[j][i]=trees[i].symmetric_difference(trees[j])
            elif matrix_type == 'euc':
                matrix[i][j]=matrix[j][i]=trees[i].euclidean_distance(trees[j])

    if invert and matrix_type == 'sym':
        matrix = max_symdiff - matrix

    if normalise:
        matrix = matrix / np.max(matrix)
        
    return matrix

def concat_distvars(outfile, filelist):
    total = len(filelist)
    current = 0
    writer = open(outfile, 'w')
    writer.write(str(total)+"\n")
    for f in filelist:
        current += 1
        contents = open(f,'r').readlines()
        writer.write( "{0} {0} {1}\n{2}".format(len(contents), current, "".join(contents) ))
    writer.close()
    return

def linkage2newick(gene_list,linkage):
    """ Generate a newick string from the cluster linkage (must supply list of taxon labels in gene_list) """
    import re
    d = {}
    num_genes = len(gene_list)
    for i in range(len(linkage)+len(gene_list)):
        if i < num_genes:
            d[i]=[gene_list[i],0]
        else:
            d[i] = [ (int(linkage[i-num_genes][0]),int(linkage[i-num_genes][1])),linkage[i-num_genes][2] ]

    s = str(i)+":0.0;"

    def quick_replace(s,d,i):
        # Replace node i with the subtree it represents
        subtree = d[i][0]
        if type(subtree)==type( (1,) ): #check tuple-ness before unpacking
            left, right = subtree
            left_brlen = d[i][1] - d[left][1] #calculate branch lengths for the unpacked subtrees from the tuple
            right_brlen = d[i][1] - d[right][1]
            out = re.sub( "(?<![\w\d\:.]){0}(?!\w\d)".format(i), "({0}:{1}, {2}:{3})".format(left,left_brlen,right,right_brlen), s ) #replace node with subtree
            return out
        else: # if tuple-ness fails then this is a leaf node, so replace with leaf name from genelist
            return re.sub( "(?<![\w\d\:.]){0}(?!\w\d)".format(i), "{0}".format(d[i][0]), s )

    # i is set to the number of the final node in the tree - this represents the whole tree contracted to a single node
    # iteratively replacing the node with the next-most contracted subtree will result in a fully expanded tree
    while i >= 0: # take advantage of i's current value at len(linkage)+len(gene_list), i.e. max
        s = quick_replace(s,d,i)
        i -= 1 # repeatedly call quick_replace function to rename node

    return s

# COLLECT FILES AND SET OPTIONS

np.set_printoptions(precision=2,linewidth=200)

print "Working directory set to {0}".format(INPUT_DIR)

taxa = dendropy.TaxonSet()
tree_files = glob.glob( "{0}/*.nwk".format(INPUT_DIR) )
if len(tree_files) == 0: 
    print "Unable to find tree files"
    sys.exit()
#info_files = glob.glob( "{0}/trees/info/*".format(INPUT_DIR) )
msa_files = glob.glob( "{0}/../../MSA/*.phy".format(INPUT_DIR) )
names = [ filename[filename.rindex("/")+1:filename.rindex(".")] for filename in tree_files]
#likelihoods = [float( re.compile( "(?<=Score of best tree ).+" ).search(open(x).read()).group() ) for x in info_files]
trees = [RAxML_object() for x in tree_files]
num_trees = len(trees)
[trees[i].read_from_path(tree_files[i],'newick',taxon_set=taxa) for i in range(num_trees)]
#for i in range(num_trees): 
#    trees[i].lnl = likelihoods[i]
num_taxa = len(taxa)
max_symdiff = 2 * (num_taxa - 3)

#Set up matrix and clustering parameters

matrix = get_distance_matrix(trees, matrix_type = matrix_type, invert = invert, normalise = normalise)
# print matrix
try: 
    Y = squareform(matrix)
    link = linkage(Y, linkage_method)
except: 
    Y = matrix
    link = linkage(Y, linkage_method)

cut = (link[-1][2])*cut
nwk = linkage2newick(names, link)


# Plot the clustering via dendrogram function - 'cut' parameter separates the trees into clusters at the given cut-point (set at 1/4 max branch length for convenience)
# T gives a list recording the number of the cluster each gene is assigned to
T = fcluster(link,cut,criterion="distance")
if make_clusters:   
    if not os.path.exists( "{0}/clusters".format(INPUT_DIR)): os.mkdir( "{0}/clusters".format(INPUT_DIR))
    if not os.path.exists( "{0}/clusters/MSA".format(INPUT_DIR)): os.mkdir( "{0}/clusters/MSA".format(INPUT_DIR))
    if not os.path.exists( "{0}/clusters/trees".format(INPUT_DIR)): os.mkdir( "{0}/clusters/trees".format(INPUT_DIR))
    clusters = {} # collect lists of sequence records for concatenation
    distvars = {} # collect separate dictionary for dv files for concatenation
    for k in range(min(T),max(T)+1):
        clusters[k] = []
        distvars[k] = []

    for i in range(len(T)):
        clusters[T[i]].append(get_phylip_file(msa_files[i]))
        distvars[T[i]].append(msa_files[i][:msa_files[i].rindex(".")]+".dv")

    for key in clusters.keys():
        seq_list = clusters[key] # list of sequence records pertaining to current cluster
        conc = concatenate_alignments(seq_list) # do concatenation of sequences
        conc.write_phylip( "{0}/clusters/MSA/cluster{1:0>2}.phy".format(INPUT_DIR,key)) # write concatenation in phylip and fasta formats
        conc.write_fasta( "{0}/clusters/MSA/cluster{1:0>2}.fas".format(INPUT_DIR,key))
        dv_list = distvars[key] # list of dv files pertaining to current cluster
        concat_distvars( "{0}/clusters/MSA/cluster{1:0>2}.dv".format(INPUT_DIR,key),dv_list) # write concatenation of dv matrices

    assignments = zip(tree_files,T)
    assignment_outfile = open( "{0}/clusters/cluster_assignments.txt".format(INPUT_DIR), "w")
    for x in sorted(assignments, key=lambda tup: int(tup[1])): 
        gene = x[0][x[0].rindex("/")+1:x[0].rindex(".")]
        assignment_outfile.write("{0}\t{1}\n".format(gene,x[1]))
    assignment_outfile.close()
    dendrogram_outfile = open( "{0}/clusters/cluster_dendrogram.nwk".format(INPUT_DIR), "w")
    dendrogram_outfile.write(nwk)
    dendrogram_outfile.close()

if show_plot:
    print nwk
    print matrix
    dendrogram( link, color_threshold=cut, leaf_font_size=10,leaf_rotation=90,leaf_label_func=lambda leaf: names[leaf]+"_"+str(T[leaf]),count_sort=True)
    title("{0} linkage of {1} matrix".format(linkage_method,matrix_type))
    axhline(cut,color='grey',ls='dashed')
    xlabel('Gene')
    ylabel('Distance')
    show()
