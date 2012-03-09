#!/usr/bin/env python
# -*- coding: utf-8 -*-

import dendropy, glob, numpy, re, sys, math, random, copy
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

def expectation_maximization(t, nbclusters=2, nbiter=3, normalize=False,\
        epsilon=0.001, monotony=False, datasetinit=True):
    """ 
    Each row of t is an observation, each column is a feature 
    'nbclusters' is the number of seeds and so of clusters
    'nbiter' is the number of iterations
    'epsilon' is the convergence bound/criterium

    Overview of the algorithm:
    -> Draw nbclusters sets of (μ, σ, P_{#cluster}) at random (Gaussian 
       Mixture) [P(Cluster=0) = P_0 = (1/n).∑_{obs} P(Cluster=0|obs)]
    -> Compute P(Cluster|obs) for each obs, this is:
    [E] P(Cluster=0|obs)^t = P(obs|Cluster=0)*P(Cluster=0)^t
    -> Recalculate the mixture parameters with the new estimate
    [M] * P(Cluster=0)^{t+1} = (1/n).∑_{obs} P(Cluster=0|obs)
        * μ^{t+1}_0 = ∑_{obs} obs.P(Cluster=0|obs) / P_0
        * σ^{t+1}_0 = ∑_{obs} P(Cluster=0|obs)(obs-μ^{t+1}_0)^2 / P_0
    -> Compute E_t=∑_{obs} log(P(obs)^t)
       Repeat Steps 2 and 3 until |E_t - E_{t-1}| < ε
    """
    def pnorm(x, m, s):
        """ 
        Compute the multivariate normal distribution with values vector x,
        mean vector m, sigma (variances/covariances) matrix s
        o in xrange(nbobs); c in xrange(nbclusters)
        Px[o,c] = pnorm(t[o,:],params[c]['mu'], params[c]['sigma'])
        """
        xmt = numpy.matrix(x-m).transpose()
        for i in xrange(len(s)):
            if s[i,i] <= sys.float_info[3]: # min float
                s[i,i] = sys.float_info[3]
        sinv = numpy.linalg.inv(s)
        xm = numpy.matrix(x-m)
        return (2.0*math.pi)**(-len(x)/2.0)*(1.0/math.sqrt(numpy.linalg.det(s)))\
                *math.exp(-0.5*(xm*sinv*xmt))

    def draw_params():
            if datasetinit:
                tmpmu = numpy.array([1.0*t[random.uniform(0,nbobs),:]],numpy.float64)
            else:
                tmpmu = numpy.array([random.uniform(min_max[f][0], min_max[f][1])\
                        for f in xrange(nbfeatures)], numpy.float64)
            return {'mu': tmpmu,\
                    'sigma': numpy.matrix(numpy.diag([(min_max[f][1]-min_max[f][0])/2.0 for f in xrange(nbfeatures)])),\
                    'proba': 1.0/nbclusters}

    nbobs = t.shape[0]
    nbfeatures = t.shape[1]
    min_max = []
    # find xranges for each features
    for f in xrange(nbfeatures):
        min_max.append((t[:,f].min(), t[:,f].max()))
    
    ### Normalization
    if normalize:
        for f in xrange(nbfeatures):
            t[:,f] -= min_max[f][0]
            t[:,f] /= (min_max[f][1]-min_max[f][0])
    min_max = []
    for f in xrange(nbfeatures):
        min_max.append((t[:,f].min(), t[:,f].max()))
    ### /Normalization

    result = {}
    quality = 0.0 # sum of the means of the distances to centroids
    random.seed()
    Pclust = numpy.ndarray([nbobs,nbclusters], numpy.float64) # P(clust|obs)
    Px = numpy.ndarray([nbobs,nbclusters], numpy.float64) # P(obs|clust)
    # iterate nbiter times searching for the best "quality" clustering
    for iteration in xrange(nbiter):
        ##############################################
        # Step 1: draw nbclusters sets of parameters #
        ##############################################
        params = [draw_params() for c in xrange(nbclusters)]
        old_log_estimate = sys.maxint         # init, not true/real
        log_estimate = sys.maxint/2 + epsilon # init, not true/real
        estimation_round = 0
        # Iterate until convergence (EM is monotone) <=> < epsilon variation
        while (abs(log_estimate - old_log_estimate) > epsilon\
                and (not monotony or log_estimate < old_log_estimate)):
            restart = False
            old_log_estimate = log_estimate
            ########################################################
            # Step 2: compute P(Cluster|obs) for each observations #
            ########################################################
            for o in xrange(nbobs):
                for c in xrange(nbclusters):
                    # Px[o,c] = P(x|c)
                    Px[o,c] = pnorm(t[o,:],\
                            params[c]['mu'], params[c]['sigma'])
            #for o in xrange(nbobs):
            #    Px[o,:] /= math.fsum(Px[o,:])
            for o in xrange(nbobs):
                for c in xrange(nbclusters):
                    # Pclust[o,c] = P(c|x)
                    Pclust[o,c] = Px[o,c]*params[c]['proba']
            #    assert math.fsum(Px[o,:]) >= 0.99 and\
            #            math.fsum(Px[o,:]) <= 1.01
            for o in xrange(nbobs):
                tmpSum = 0.0
                for c in xrange(nbclusters):
                    tmpSum += params[c]['proba']*Px[o,c]
                Pclust[o,:] /= tmpSum
                #assert math.fsum(Pclust[:,c]) >= 0.99 and\
                #        math.fsum(Pclust[:,c]) <= 1.01
            ###########################################################
            # Step 3: update the parameters (sets {mu, sigma, proba}) #
            ###########################################################
            print "iter:", iteration, " estimation#:", estimation_round,\
                    " params:", params
            for c in xrange(nbclusters):
                tmpSum = math.fsum(Pclust[:,c])
                params[c]['proba'] = tmpSum/nbobs
                if params[c]['proba'] <= 1.0/nbobs:           # restart if all
                    restart = True                             # converges to
                    print "Restarting, p:",params[c]['proba'] # one cluster
                    break
                m = numpy.zeros(nbfeatures, numpy.float64)
                for o in xrange(nbobs):
                    m += t[o,:]*Pclust[o,c]
                params[c]['mu'] = m/tmpSum
                s = numpy.matrix(numpy.diag(numpy.zeros(nbfeatures, numpy.float64)))
                for o in xrange(nbobs):
                    s += Pclust[o,c]*(numpy.matrix(t[o,:]-params[c]['mu']).transpose()*\
                            numpy.matrix(t[o,:]-params[c]['mu']))
                    #print ">>>> ", t[o,:]-params[c]['mu']
                    #diag = Pclust[o,c]*((t[o,:]-params[c]['mu'])*\
                    #        (t[o,:]-params[c]['mu']).transpose())
                    #print ">>> ", diag
                    #for i in xrange(len(s)) :
                    #    s[i,i] += diag[i]
                params[c]['sigma'] = s/tmpSum
                print "------------------"
                print params[c]['sigma']

            ### Test bound conditions and restart consequently if needed
            if not restart:
                restart = True
                for c in xrange(1,nbclusters):
                    if not numpy.allclose(params[c]['mu'], params[c-1]['mu'])\
                    or not numpy.allclose(params[c]['sigma'], params[c-1]['sigma']):
                        restart = False
                        break
            if restart:                # restart if all converges to only
                old_log_estimate = sys.maxint          # init, not true/real
                log_estimate = sys.maxint/2 + epsilon # init, not true/real
                params = [draw_params() for c in xrange(nbclusters)]
                continue
            ### /Test bound conditions and restart

            ####################################
            # Step 4: compute the log estimate #
            ####################################
            log_estimate = math.fsum([math.log(math.fsum(\
                    [Px[o,c]*params[c]['proba'] for c in xrange(nbclusters)]))\
                    for o in xrange(nbobs)])
            print "(EM) old and new log estimate: ",\
                    old_log_estimate, log_estimate
            estimation_round += 1

        # Pick/save the best clustering as the final result
        quality = -log_estimate
        if not quality in result or quality > result['quality']:
            result['quality'] = quality
            result['params'] = copy.deepcopy(params)
            result['clusters'] = [[o for o in xrange(nbobs)\
                    if Px[o,c] == max(Px[o,:])]\
                    for c in xrange(nbclusters)]
    return result

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
s = expectation_maximization(matrix,4,6,normalize=False,monotony=False,datasetinit=True)
print s
