#!/usr/bin/env 

def handleArgs(argv,help="no arguments given\n"):
    import sys # N.B. implicitly sys should already be loaded in order to use this function
    """ argv = sys.argv, i.e a list of commandline arguments, 
        with -a -b type flags
        It returns a dictionary in which flags are keys and their following arguments
        are values. If a flag is not followed by a value the dictionary stores 'None'
        Flags can be given in any order, and repetition overwrites
    """
    if len(argv) < 2: 
        print help
        sys.exit()
            
    argv = argv[1:]
    flagDic = {}
    for i in range(len(argv)):
        if argv[i].startswith("-"):
            try: 
                if not argv[i+1].startswith("-"): 
                    flagDic[argv[i]] = argv[i+1]
                else: flagDic[argv[i]] = None
            except IndexError: flagDic[argv[i]] = None
    return flagDic

def pam2sps(tree_file, conversion, outfile=None):
    """ Take in a newick file, change branch lengths by factor 100, write output file
        (i.e. convert between PAM units and substitutions per site) """
    import re
    reg_ex = re.compile('(?<=:)[0-9.]+')
    convert_pam_to_sps = lambda a: str(0.01*float(a.group()))
    convert_sps_to_pam = lambda b: str(100*float(b.group()))

    input_string = open(tree_file).read()
    if conversion == 'pam2sps':
        output_string = reg_ex.sub(convert_pam_to_sps,input_string)
    elif conversion == 'sps2pam':
        output_string = reg_ex.sub(convert_sps_to_pam,input_string)
    else: output_string = reg_ex.sub(convert_pam_to_sps,input_string)

    if outfile:
        writer = open(outfile,'w')
        writer.write(output_string)
        writer.close()
    else: print output_string

def get_alignments(input_dir):
    # Returns a list of fasta files in a given directory (provided the file extension is .fa or .fas)
    # Doesn't allow for a mix of .fa and .fas
    import glob
    fasta_files = glob.glob( "{0}/*.fa".format(input_dir) )
    if len(fasta_files) == 0:
        fasta_files = glob.glob( "{0}/*.fas".format(input_dir) )
    fasta_files = sorted(fasta_files)
    return fasta_files

def get_gene_trees(input_dir):
    # Returns a list of Inference_Result objects from files in input_dir
    import glob, inference_functions
    result_files = glob.glob( "{0}/*.trobj".format(input_dir) ) # using file ext 'trobj' for now
    for result_file in result_files:
        result = inference_functions.Inference_Result()
        result.read_from_file(result_file)
        yield result

def populate_phylip_from_fasta(fastafile):
    import sequence_record
    phylip_name = fastafile[:fastafile.rindex(".")] + ".phy"
    sequence_record.get_fasta_file(fastafile).write_phylip(phylip_name)

def populate_dv_from_fasta(fastafile, datatype, helper="./library_darwin_code/TC_wrapper.drw"):
    import os
    prefix = fastafile[:fastafile.rindex(".")]
    command = 'echo "fil := ReadFastaWithNames(\'{0}\'); seqtype := \'{1}\' ; ReadProgram(\'{2}\');" | darwin > /dev/null'.format(fastafile, datatype, helper)
    os.system( command )
    os.rename( "temp_distvar.txt", prefix+"_dv.txt")
    os.rename( "temp_map.txt", prefix+"_map.txt")
    os.rename( "temp_labels.txt", prefix+"_labels.txt")
    os.rename( "temp_tree.nwk", prefix+"_guidetree.nwk")

############
# Concatenating TreeCollection input files:
# Two helper functions and the main one

def _parse_dv(distvar_file):
    raw = [line.rstrip() for line in open(distvar_file).readlines()]
    while '' in raw: raw.remove('')
    num_matrices = int(raw.pop(0))
    matrices = []
    read_records = 0
    while read_records < num_matrices:
        dims = raw.pop(0).split()
        size = int(dims[0])
        index = dims[-1]
        c = 0
        matrix = ''
        while c < size:
            matrix += raw.pop(0) + "\n"
            c += 1
        matrices.append(matrix)
        read_records += 1
    return matrices

def _parse_maps(map_file):
    raw = open(map_file).readlines()
    while '' in raw: raw.remove('')
    header = raw.pop(0)
    maps = ''.join(raw)
    return maps
    
def concatenate_dvs(distvar_files, map_files, labels, tree, outprefix=None):
    input_files = zip(distvar_files, map_files)
    dvstring    = ''
    mapstring   = ''
    index       = 0
    for d,m in input_files:
        dvmatrices = _parse_dv(d)
        for dv in dvmatrices:
            index     += 1
            dim        = "{0} {0} {1}\n".format( len(dv.rstrip().split("\n") ), index )
            dvstring  += (dim+dv)
            mapstring += _parse_maps(m)
    dvstring     = "{0}\n{1}".format(index, dvstring)
    mapstring    = "{0} {1}\n{2}".format(index, open(labels).readline().rstrip(),mapstring)
    labelsstring = open(labels).read()
    treestring   = open(tree).read()
    if outprefix:
        dv_out     = outprefix + "_dv.txt"
        map_out    = outprefix + "_map.txt"
        labels_out = outprefix + "_labels.txt"
        tree_out   = outprefix + "_guidetree.nwk"
        open(dv_out,'w').write(dvstring)
        open(map_out,'w').write(mapstring)
        open(labels_out,'w').write(labelsstring)
        open(tree_out,'w').write(treestring)

############

def get_distance_matrix(trees, matrix_type="sym", invert=False, normalise=False):
    import numpy as np
    import dendropy as dpy
    """ all pairwise distances between trees in list: 'trees' """
    np.set_printoptions(precision=2,linewidth=200)
    num_trees = len(trees)
    matrix = np.zeros( (num_trees,num_trees),dtype='float' )
    taxa = dpy.TaxonSet()
    dpytrees = [dpy.Tree() for tree in trees]
    for i in range(num_trees):
        dpytrees[i].read_from_string(trees[i].tree,'newick',taxon_set=taxa)
    for i in range(num_trees):
        for j in range(i+1,num_trees):
            if matrix_type == 'rf':
                matrix[i][j]=matrix[j][i]=dpytrees[i].robinson_foulds_distance(dpytrees[j])
            elif matrix_type == 'sym':
                matrix[i][j]=matrix[j][i]=dpytrees[i].symmetric_difference(dpytrees[j])
            elif matrix_type == 'euc':
                matrix[i][j]=matrix[j][i]=dpytrees[i].euclidean_distance(dpytrees[j])

    if invert and matrix_type == 'sym':
        max_symdiff = 2 * (len(taxa) - 3)
        matrix = max_symdiff - matrix

    if normalise:
        matrix = matrix / np.max(matrix)
  
    return matrix

def get_linkage(matrix, linkage_method="ward"):
    """ Linkage methods are: single, complete, average, ward """
    from scipy.cluster.hierarchy import linkage
    from hcluster import squareform
    try: 
        Y = squareform(matrix)
        link = linkage(Y, linkage_method)
    except: 
        Y = matrix
        link = linkage(Y, linkage_method)
    return link

def cluster_linkage(link, threshold, criterion="maxclust"):
    """ Returns list of cluster assignments from linkage
        Criteria: 'maxclust' - set threshold as (maximum) number of groups to cluster into 
                  'distance' - set threshold as cutpoint at which to separate dendrogram
                               into clusters (threshold in range float(0,1)) """
    from scipy.cluster.hierarchy import fcluster, maxinconsts
    if criterion == "distance":
        threshold = (link[-1][2])*threshold
    T = fcluster(link, threshold, criterion=criterion)
    return T

def assign_to_clusters(msa_files, T, output_dir=None):
    from sequence_record import *
    import os
    clusters    = {} # collect lists of sequence records for concatenation
    distvars    = {} # collect separate dictionary for dv files for concatenation
    maps        = {}
    assignments = {}
    for k in range(min(T),max(T)+1):
        clusters[k]    = []
        distvars[k]    = []
        maps[k]        = []
        assignments[k] = []
    for i in range(len(T)):
        clusters[T[i]].append(get_phylip_file(msa_files[i]))
        distvars[T[i]].append(msa_files[i][:msa_files[i].rindex(".")]+"_dv.txt")
        maps[T[i]].append(msa_files[i][:msa_files[i].rindex(".")]+"_map.txt")
        assignments[T[i]].append(msa_files[i][msa_files[i].rindex("/")+1:msa_files[i].rindex(".")])

    for key in clusters:
        seq_list   = clusters[key] # list of sequence records pertaining to current cluster
        dv_list    = distvars[key] # list of dv files pertaining to current cluster
        map_list   = maps[key]
        labels     = dv_list[0][:dv_list[0].rindex("_")] + "_labels.txt"
        guidetree = dv_list[0][:dv_list[0].rindex("_")] + "_guidetree.nwk"
        conc       = concatenate_alignments(seq_list) # do concatenation of sequences
        if output_dir:
            if not os.path.isdir(output_dir): os.mkdir(output_dir)
            conc.write_phylip("{0}/cluster{1:0>2}.phy".format(output_dir, key))
            conc.write_fasta("{0}/cluster{1:0>2}.fas".format(output_dir, key))
            conc_dvs   = concatenate_dvs(dv_list, map_list, labels, guidetree, "{0}/cluster{1:0>2}".format(output_dir,key))

    return assignments
        
##################################################################################################
