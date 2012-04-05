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

def check_args_filepath( arg, argDic, alternative ):
    import os
    if arg in argDic: 
        if argDic[arg]:
            filepath = argDic[arg]
            while filepath.endswith("/"): filepath = filepath[:-1] # Trim off trailing "/"s
        else: filepath = alternative
    else: filepath = alternative
    if filepath:  
        if not os.path.isdir( filepath):
            os.mkdir( filepath )
    return filepath

def check_args_bool( arg, argDic ):
    if arg in argDic: 
        return True
    else:
        return False

def check_args_value( arg, argDic, alternative ):
    if arg in argDic: 
        if argDic[arg]:
            value = argDic[arg]
        else: value = alternative
    else: value = alternative
    return value

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

def get_alignments(input_dir, format='fasta'):
    # Returns a list of fasta files in a given directory (provided the file extension is .fa or .fas)
    # Doesn't allow for a mix of .fa and .fas
    import glob
    if format=='fasta':
        fasta_files = glob.glob( "{0}/*.fa".format(input_dir) )
        if len(fasta_files) == 0:
            fasta_files = glob.glob( "{0}/*.fas".format(input_dir) )
        fasta_files = sorted(fasta_files)
        return fasta_files
    elif format=='dv':
        dv_files = glob.glob( "{0}/*dv.txt".format(input_dir))
        return dv_files
    else: return []

def get_gene_trees(input_dir):
    # Returns a list of Inference_Result objects from files in input_dir
    import glob, inference_functions
    result_files = glob.glob( "{0}/*.tree".format(input_dir) ) # using file ext '.tree' for Inference_Result objects (also standard newick compatible)
    for result_file in result_files:
        result = inference_functions.Inference_Result()
        result.read_from_file(result_file)
        yield result

def populate_phylip_from_fasta(fastafile):
    import sequence_record
    phylip_name = fastafile[:fastafile.rindex(".")] + ".phy"
    sequence_record.get_fasta_file(fastafile).write_phylip(phylip_name)

def sanitise_fasta(fasta):
    import sequence_record as sr
    s = sr.get_fasta_file(fasta)
    s.sort_by_name()
    l = []
    for h in s.headers:
        while h.startswith(' '): h = h[1:]
        h = h.replace(' ','_')
        l.append(h)
    s.headers=l
    s.write_fasta(fasta)

def populate_dv_from_fasta(fastafile, datatype, helper="./library_darwin_code/TC_wrapper.drw", fpath="./"):
    import os, re
    from subprocess import Popen, PIPE
    if not fpath.endswith("/"): fpath += "/"
    if not os.path.isdir(fpath): os.mkdir(fpath)
    prefix = fastafile[:fastafile.rindex(".")]
    command = 'echo "fil := ReadFastaWithNames(\'{0}\'); seqtype := \'{1}\'; fpath := \'{2}\'; ReadProgram(\'{3}\');" | darwin'.format(fastafile, datatype, fpath, helper)
    process = Popen( command, shell=True, stdout=PIPE, stderr=PIPE )
    stdout,stderr = process.communicate()
    score = float(re.search("(?<=Score: )[0-9.]+",stdout).group())
    tree = re.search("(?<=Tree: ).+(?=\n)",stdout).group()
    labels = re.search("(?<=Labels: ).+(?=\n)",stdout).group()
    os.rename( "temp_distvar.txt", prefix+"_dv.txt")
    os.rename( "temp_map.txt", prefix+"_map.txt")
    os.rename( "temp_labels.txt", prefix+"_labels.txt")
    os.rename( "temp_tree.nwk", prefix+"_guidetree.nwk")
    return score, tree, labels

############
# Concatenating TreeCollection input files:
# Two helper functions and the main one
    
def concatenate_dvs(distvar_files, map_files, labels, tree, outprefix=None):
    
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
    dpytrees = [dpy.Tree.get_from_string(tree.tree,'newick',taxon_set=taxa) for tree in trees]
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
    except ValueError: 
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
        link_size = len(link)
        if threshold <= 1:
            br_top = link[link_size-threshold][2]
        else:
            br_top = link[link_size-threshold+1][2]
        if threshold  >= len(link):
            br_bottom = 0
        else: 
            br_bottom = link[link_size-threshold][2]
        threshold = 0.5 * (br_top + br_bottom)
    T = fcluster(link, threshold, criterion=criterion)

    return T

def assign_to_clusters(msa_files, T,output_dir=None):
    import sequence_record as sr
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
        clusters[T[i]].append( sr.get_phylip_file(msa_files[i]) )
        assignments[T[i]].append(msa_files[i][msa_files[i].rindex("/")+1:msa_files[i].rindex(".")])
        distvars[T[i]].append(msa_files[i][:msa_files[i].rindex(".")]+"_dv.txt")
        maps[T[i]].append(msa_files[i][:msa_files[i].rindex(".")]+"_map.txt")

    for key in clusters:
        seq_list   = clusters[key] # list of sequence records pertaining to current cluster
        conc       = sr.concatenate_alignments(seq_list) # do concatenation of sequences
        dv_list    = distvars[key] # list of dv files pertaining to current cluster
        map_list   = maps[key]
        labels     = dv_list[0][:dv_list[0].rindex("_")] + "_labels.txt"
        guidetree = dv_list[0][:dv_list[0].rindex("_")] + "_guidetree.nwk"
    
        if output_dir:
            if not os.path.isdir(output_dir): os.mkdir(output_dir)
            conc.write_phylip("{0}/cluster{1:0>2}.phy".format(output_dir, key))
            conc.write_fasta("{0}/cluster{1:0>2}.fas".format(output_dir, key))
            conc_dvs = concatenate_dvs(dv_list, map_list, labels, guidetree, "{0}/cluster{1:0>2}".format(output_dir,key))

    return assignments

def assign_to_clusters_optimiser(msa_files, T, output_dir=None):
    import sequence_record as sr
    import os
    distvars    = {} # collect separate dictionary for dv files for concatenation
    maps        = {}
    assignments = {}
    for k in range(min(T),max(T)+1):
        distvars[k]    = []
        maps[k]        = []
        assignments[k] = []
    for i in range(len(T)):
        assignments[T[i]].append(msa_files[i][msa_files[i].rindex("/")+1:msa_files[i].rindex(".")])
        distvars[T[i]].append(msa_files[i][:msa_files[i].rindex(".")]+"_dv.txt")
        maps[T[i]].append(msa_files[i][:msa_files[i].rindex(".")]+"_map.txt")

    for key in assignments:
        dv_list    = distvars[key] # list of dv files pertaining to current cluster
        map_list   = maps[key]
        try: 
            labels     = dv_list[0][:dv_list[0].rindex("_")] + "_labels.txt"
            guidetree = dv_list[0][:dv_list[0].rindex("_")] + "_guidetree.nwk"
        except IndexError: continue
    
        if output_dir:
            if not os.path.isdir(output_dir): os.mkdir(output_dir)
            conc_dvs = concatenate_dvs(dv_list, map_list, labels, guidetree, "{0}/cluster{1:0>2}".format(output_dir,key))

    return assignments

def calc_distinct_groups(matrix):
    nclusters = len(matrix)
    indices = range(nclusters)
    for i in range(nclusters):
        if i in indices:
            for j in range(i+1, nclusters):
                if matrix[i][j] == 0:
                    indices.remove(j)
    return len(indices)

def showplot(matrix, T, link, names, nclasses=None):
    import scipy.cluster.hierarchy as hchy
    import matplotlib.pyplot as plt
    if nclasses: 
        link_size = len(link)
        if nclasses <= 1:
            br_top = link[link_size-nclasses][2]+0.01
        else:
            br_top    = link[link_size-nclasses+1][2]
        if nclasses  >= len(link):
            br_bottom = 0
        else: 
            br_bottom = link[link_size-nclasses][2]
        cut = 0.5 * (br_top + br_bottom)
    else: cut = (link[-1][2])*0.25
    hchy.dendrogram( link, color_threshold=cut, leaf_font_size=10,leaf_rotation=90,leaf_label_func=lambda leaf: names[leaf]+"_"+str(T[leaf]),count_sort=True)
    plt.title("Dendrogram")
    plt.axhline(cut,color='grey',ls='dashed')
    plt.xlabel('Gene')
    plt.ylabel('Distance')
    plt.show()
    plt.clf()

def compute_constrained_score( dist_mat, var_mat, labels, cluster_tree, helper="./library_darwin_code/conscore.drw"):
    from subprocess import PIPE, Popen
    import os, re
    import numpy as np
    np.savetxt( "tmp_dv.txt", dist_mat, fmt='%.2f')
    np.savetxt( "tmp_var.txt", var_mat, fmt='%.4f')
    cmd = ' echo " dm_raw_text := ReadRawFile(\'tmp_dv.txt\'): labels := [{0}]: vm_raw_text := ReadRawFile(\'tmp_var.txt\'): tree := ParseNewickTree(\'{1}\'): ReadProgram(\'{2}\'): " | darwin '.format( labels, cluster_tree, helper )
    process = Popen( cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    #print stdout
    scores = re.search("(?<=fitting index )([0-9.]+)\n([0-9.]+)\n([0-9.]+)",stdout)
    os.remove("tmp_dv.txt")
    os.remove("tmp_var.txt")
    return [float(x) for x in scores.groups()] # Dimensionless fitting index,  Sum of Squares,  MST_Qual (sum of squares / degrees of freedom)

def compute_weighted_robinson_foulds( gene_tree, cluster_tree ):
    import dendropy as dpy
    g = dpy.Tree().get_from_string(gene_tree, 'newick')
    c = dpy.Tree().get_from_string(cluster_tree, 'newick')
    return g.robinson_foulds_distance(c)

def compute_robinson_foulds( gene_tree, cluster_tree ):
    import dendropy as dpy
    g = dpy.Tree().get_from_string(gene_tree, 'newick')
    c = dpy.Tree().get_from_string(cluster_tree, 'newick')
    return g.symmetric_difference(c)

def separate_dv_matrix( dv_file ):
    import numpy as np
    np.set_printoptions(suppress=True, precision=2)
    dv = open( dv_file ).readlines()[1:]
    length = int(dv.pop(0).split()[0])
    a = np.zeros( (length, length) )
    i = 0
    while dv:
        row = [float(n) for n in dv.pop(0).rstrip().split()]
        a[i] = row
        i += 1
    d = np.triu(a)
    v = np.tril(a)

    for i in range(length):
        v[i] = v.T[i]
        d.T[i] = d[i]

    return d,v

def draw_rand():
    import random
    return random.random()

def probs(l):
    import math
    l2 = [math.exp(-x) for x in l]
    s = sum(l2)
    prev = 0
    p = []
    for v in l2:
        v = v/s + prev
        p.append(v)
        prev = v
    return p

def probs2(l):
    l2 = [1/x for x in l]
    s = sum(l2)
    prev = 0
    p = []
    for v in l2:
        v = v/s + prev
        p.append(v)
        prev = v
    return p

def choice(move_probs):
    import random
    n = random.random()
    for k in range(len(move_probs)):
        if n > move_probs[k]: continue
        else: 
            move_to = k+1
            break
    return move_to

def get_cluster_trees( cluster_tree_dir, quiet=False ):
    import glob
    from inference_functions import run_treecollection
    cluster_files = glob.glob( "{0}/*_dv.txt".format(cluster_tree_dir) )
    cluster_names = [ x[x.rindex("/")+1:x.rindex("_")] for x in cluster_files]
    cluster_trees = []
    for cluster in cluster_files:
        prefix = cluster[:cluster.rindex("_")]
        dv = prefix+"_dv.txt"       
        map_file = prefix+"_map.txt"
        labels = prefix+"_labels.txt"
        guide_tree = prefix+"_guidetree.nwk"
        name = prefix[prefix.rindex("/")+1:]
        if not quiet: print "Running TreeCollection on {0}".format(dv)
        tree = run_treecollection(dv, map_file, labels, guide_tree, name)
        cluster_trees.append( tree )
    return cluster_trees

def optimise_sample(clustering, gene_trees, sample_size, max_reassignments, INPUT_DIR, CLUSTER_DIR, fasta_files, best_score, greedy=False):
    import copy, random
    new_clustering = copy.copy(clustering)
    dic = {}
    sample = random.sample(range(len(gene_trees)), sample_size)
    for i in sample:
        test_tree = gene_trees[i]
        dv_file  = "{0}/{1}_dv.txt".format(INPUT_DIR,test_tree.name)
        labels = open("{0}/{1}_labels.txt".format(INPUT_DIR,test_tree.name)).readlines()[1].rstrip()
        labels = ', '.join(["'{0}'".format(x) for x in labels.split()])
        d,v = separate_dv_matrix( dv_file )
        belongs_to = clustering[i]
        scores = []
        for j in range(len(cluster_trees)):
            topology = cluster_trees[j].tree
            score = compute_constrained_score( d, v, labels, topology )[1]
            scores.append( score )
        if belongs_to != choice(probs(scores)):
            dic[i] = (choice(probs(scores)))   
    for k in random.sample(dic.keys(),min(max_reassignments,len(dic))):
        new_clustering[k] = dic[k]
    assignments = assign_to_clusters_optimiser(fasta_files, new_clustering, CLUSTER_DIR)
    new_cluster_trees = get_cluster_trees(CLUSTER_DIR, quiet=True)
    new_score = sum([float(tr.score) for tr in new_cluster_trees])
    
    # Acceptance decision
    if greedy:
        if new_score < best_score: # greedy
            return new_clustering, new_score
        else:
            os.system( "rm {0}/*".format(CLUSTER_DIR) )
            assignments = assign_to_clusters_optimiser(fasta_files, clustering, CLUSTER_DIR)
            return clustering, best_score
        
    else:
        if random.random() > (new_score / (new_score + best_score)): # stochastic
            return new_clustering, new_score
        else:
            os.system( "rm {0}/*".format(CLUSTER_DIR) )
            assignments = assign_to_clusters_optimiser(fasta_files, clustering, CLUSTER_DIR)
            return clustering, best_score

def optimise_sample_rf(clustering, gene_trees, cluster_trees, sample_size, max_reassignments, INPUT_DIR, CLUSTER_DIR, fasta_files, best_score, greedy=False):
    import copy, os, random
    new_clustering = copy.copy(clustering)
    dic = {}
    sample = random.sample(range(len(gene_trees)), sample_size)
    for i in sample:
        test_tree = gene_trees[i]
        belongs_to = clustering[i]
        scores = []
        for j in range(len(cluster_trees)):
            topology = cluster_trees[j].tree
            score = compute_robinson_foulds(test_tree.tree, topology)
            scores.append( score )
        moves_to = choice(probs(scores))
        #print scores, probs(scores)
        if belongs_to != moves_to:
            dic[i] = moves_to
    for k in random.sample(dic.keys(),min(max_reassignments,len(dic))):
        new_clustering[k] = dic[k]
    assignments = assign_to_clusters_optimiser(fasta_files, new_clustering, CLUSTER_DIR)
    new_cluster_trees = get_cluster_trees(CLUSTER_DIR, quiet=True)
    new_score = sum([float(tr.score) for tr in new_cluster_trees])
    
    # Acceptance decision
    if greedy:
        if new_score < best_score: # greedy
            return new_clustering, new_score
        else:
            os.system( "rm {0}/*".format(CLUSTER_DIR) )
            assignments = assign_to_clusters_optimiser(fasta_files, clustering, CLUSTER_DIR)
            return clustering, best_score
        
    else:
        if random.random() > (new_score / (new_score + best_score)): # stochastic
            return new_clustering, new_score
        else:
            os.system( "rm {0}/*".format(CLUSTER_DIR) )
            assignments = assign_to_clusters_optimiser(fasta_files, clustering, CLUSTER_DIR)
            return clustering, best_score

##################################################################################################
