#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import os
import dendropy as dpy
import GeoMeTreeHack
from result import Result
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from Bio.Cluster import kmedoids
import matplotlib.pyplot as plt
from matplotlib import cm as CM
from math import log
from copy import copy, deepcopy
from sklearn.cluster import KMeans, spectral_clustering, \
    AffinityPropagation
from mpl_toolkits.mplot3d import Axes3D
from collections import defaultdict


class Clustering(object):

    """
    Class to store clustering information
    """

    def __init__(self):
        """
        self.partitions uses a compound key to retrieve partitions
        key = tuple of (distance_metric, linkage_method, num_classes)
        """

        self.distance_matrices = {}
        self.partitions = {}
        self.plotting_info = {}
        self.clusters = {}
        # self._affinity_propagation = {
        #     'euc': defaultdict(list),
        #     'geodesic': defaultdict(list),
        #     'sym': defaultdict(list),
        #     'rf': defaultdict(list),
        #     }

    def __str__(self):
        s = ''
        num_partitions = len(self.partitions)
        num_clusters = len(self.clusters)
        if num_clusters == 1:
            s += '{0} partition calculated:\n'.format(num_partitions)
        else:
            s += '{0} partitions calculated:\n'.format(num_partitions)
        for p in self.partitions:
            s += ' '.join(str(x) for x in p) + '\n'
        return s

    def get_distance_matrix(
        self,
        trees,
        metric,
        invert=False,
        normalise=False,
        tmpdir='/tmp',
        gtp_path='./class_files',
        ):
        """
        Generates pairwise distance matrix between trees
        Uses one of the following distance metrics:
        Robinson-Foulds distance - topology only (='sym')
        Robinson-Foulds distance - branch lengths (='rf')
        Euclidean distance - Felsenstein's branch lengths
            distance (='euc')
        Geodesic distance - branch lengths (='geodesic')
        """

        num_trees = len(trees)
        matrix = np.zeros((num_trees, num_trees), dtype='float')
        if metric == 'geometree':
            geotrees = [tree.newick for tree in trees]
            for i in range(num_trees):
                for j in range(i + 1, num_trees):
                    matrix[i][j] = matrix[j][i] = \
                        GeoMeTreeHack.main(geotrees[i], geotrees[j])
        elif metric == 'geo':
            rooted = all([tree.rooted for tree in trees])
            with open('{0}/geotrees.nwk'.format(tmpdir),'w') as file:            
                file.write('\n'.join([tree.newick.rstrip() for tree in trees]))
            if rooted:
                print 'All trees are rooted'
                os.system('java -jar {1}/gtp.jar -o {0}/output.txt {0}/geotrees.nwk'.format(tmpdir, gtp_path))
            else:
                print 'Not all trees are rooted'
                os.system('java -jar {1}/gtp.jar -u -o {0}/output.txt {0}/geotrees.nwk'.format(tmpdir, gtp_path))
            with open('{0}/output.txt'.format(tmpdir)) as file:
                for line in file:
                    line = line.rstrip()
                    if line:
                        i,j,value=line.split()
                        i = int(i)
                        j = int(j)
                        value=float(value)
                        matrix[i,j] = matrix[j,i] = value
            os.remove('{0}/output.txt'.format(tmpdir))
            os.remove('{0}/geotrees.nwk'.format(tmpdir))
        else:
            taxa = dpy.TaxonSet()
            dpytrees = [dpy.Tree.get_from_string(tree.newick, 'newick',
                        taxon_set=taxa) for tree in trees]
            for i in range(num_trees):
                for j in range(i + 1, num_trees):
                    if metric == 'rf':
                        matrix[i][j] = matrix[j][i] = \
                            dpytrees[i].robinson_foulds_distance(dpytrees[j])
                    elif metric == 'sym':
                        dist = dpytrees[i].symmetric_difference(dpytrees[j])
                        matrix[i][j] = matrix[j][i] = dist
                    elif metric == 'euc':
                        matrix[i][j] = matrix[j][i] = \
                            dpytrees[i].euclidean_distance(dpytrees[j])

            if invert and metric == 'sym':
                max_symdiff = 2 * (len(taxa) - 3)
                matrix = max_symdiff - matrix

        if normalise:
            matrix = matrix / np.max(matrix)

        return matrix

    def put_distance_matrix(
        self,
        trees,
        metric,
        invert=False,
        normalise=False,
        tmpdir='/tmp',
        gtp_path='./class_files'
        ):
    	matrix = self.get_distance_matrix(trees, metric, invert=invert,
                normalise=normalise, tmpdir=tmpdir, gtp_path=gtp_path)
        self.distance_matrices[metric] = matrix

    def order(self, l):
        """
        The clustering returned by the hcluster module gives 
        group membership without regard for numerical order 
        This function preserves the group membership, but sorts 
        the labelling into numerical order
        
        *Faster version without recursion*
        """

        list_length = len(l)

        d = defaultdict(list)
        for (i, element) in enumerate(l):
            d[element].append(i)

        l2 = [None]*list_length

        for (name, index_list) in enumerate(sorted(d.values(), key=min), start=1):
            for index in index_list:            
                l2[index] = name

        return l2

    def put_partition(
        self,
        metric,
        linkage_method,
        nclasses,
        names,
        criterion='distance',
        prune=True,
        ):
        """ Returns list of cluster assignments from linkage
            Criteria: 'maxclust' -  set threshold as (maximum) number 
                                    of groups to cluster into 
                      'distance' -  set threshold as cutpoint at which 
                                    to separate dendrogram into clusters 
                                    (threshold in range float(0,1)) """

        dm = self.distance_matrices[metric]
        compound_key = (metric, linkage_method, nclasses)
        if nclasses == 1:
            T = [1] * len(names)
            self.partitions[compound_key] = T
            return

        if linkage_method == 'kmedoids':
            p = []
            for i in range(100):
                p.append(kmedoids(dm, nclusters=nclasses, npass=100))
            T = sorted(p, key=lambda x: x[1])[0][0]
        elif linkage_method == 'MDS':

            dbc = self.get_double_centre(dm)
            (eigvals, eigvecs, cve) = self.get_eigen(dbc, standardize=True)
            coords = self.get_coords_by_cutoff(eigvals, eigvecs, cve,
                    95, normalise=False)
            est = KMeans(n_clusters=nclasses)
            est.fit(coords)
            T = est.labels_
        elif linkage_method == 'spectral':

            laplacian = self.spectral(dm, prune=prune)
            (eigvals, eigvecs, cve) = self.get_eigen(laplacian, standardize=False)
            coords = self.get_coords_by_dimension(eigvals, eigvecs,
                    cve, nclasses, normalise=True)[0]
            est = KMeans(n_clusters=nclasses)
            est.fit(coords)
            T = est.labels_
        elif linkage_method == 'affinity':

            T = self.affinity_propagation(dm, metric, nclasses)
        else:
            if metric == 'sym':
                size = len(dm)
                new = np.zeros( (size,size) )
                for i in range(size):
                    for j in range(i+1,size):
                        eps = np.random.normal(0,0.001)
                        if dm[i,j] + eps > 0:
                            new[i,j]=new[j,i]=dm[i,j]+eps
                        else:
                            new[i,j]=new[j,i]=dm[i,j]-eps
                dm = new
            linkmat = linkage(dm, linkage_method)

            if criterion == 'distance':
                linkmat_size = len(linkmat)
                if nclasses <= 1:
                    br_top = linkmat[linkmat_size - nclasses][2]
                else:
                    br_top = linkmat[linkmat_size - nclasses + 1][2]
                if nclasses >= len(linkmat):
                    br_bottom = 0
                else:
                    br_bottom = linkmat[linkmat_size - nclasses][2]
                threshold = 0.5 * (br_top + br_bottom)
            T = fcluster(linkmat, threshold, criterion=criterion)
            self.plotting_info[compound_key] = (linkmat, names,
                    threshold)

        T = self.order(T)  # puts cluster labels in ascending order
        self.partitions[compound_key] = T

    def plot_embedding(
        self,
        metric,
        linkage,
        nclasses,
        dimensions=3,
        embed='MDS',
        standardize=False,
        normalise=True
        ):

        if not dimensions in [2, 3]:
            print '2D or 3D only'
            return

        dm = self.distance_matrices[metric]
        partition = self.partitions[(metric, linkage, nclasses)]

        if embed == 'MDS':
            dbc = self.get_double_centre(dm)
            (eigvals, eigvecs, cve) = self.get_eigen(dbc, standardize=standardize)
            (coords, varexp) = self.get_coords_by_dimension(eigvals,
                    eigvecs, cve, dimensions, normalise=normalise)
        elif embed == 'spectral':

            laplacian = self.spectral(dm)
            (eigvals, eigvecs, cve) = self.get_eigen(laplacian, standardize=standardize)
            (coords, varexp) = self.get_coords_by_dimension(eigvals,
                    eigvecs, cve, dimensions, normalise=normalise)
        else:

            print 'Embedding must be one of MDS or spectral (default=MDS)'
            return

        colors = 'bgrcmyk'
        fig = plt.figure()
        if dimensions == 3:
            ax = fig.add_subplot(111, projection='3d')
        else:
            ax = fig.add_subplot(111)
        for i in range(len(partition)):
            ax.scatter(color=colors[partition[i] % len(colors)],
                       *coords[i])
        return fig

    def plot_dendrogram(self, compound_key):
        """
        Extracts data from clustering to plot dendrogram
        """

        partition = self.partitions[compound_key]
        (linkmat, names, threshold) = self.plotting_info[compound_key]
        fig = plt.figure(figsize=(11.7, 8.3))
        dendrogram(
            linkmat,
            color_threshold=threshold,
            leaf_font_size=8,
            leaf_rotation=90,
            leaf_label_func=lambda leaf: names[leaf] + '_' \
                + str(partition[leaf]),
            count_sort=True,
            )
        plt.suptitle('Dendrogram', fontsize=16)
        plt.title('Distance metric: {0}    Linkage method: {1}    Number of classes: {2}'.format(compound_key[0],
                  compound_key[1], compound_key[2]), fontsize=12)
        plt.axhline(threshold, color='grey', ls='dashed')
        plt.xlabel('Gene')
        plt.ylabel('Distance')
        return fig

    def plot_heatmap(self, metric, gridsize=None):

        dm = self.distance_matrices[metric]
        length = len(dm)
        X, Y = np.meshgrid( np.arange(length),
                            np.arange(length))
        x = X.ravel()
        y = Y.ravel()
        z = dm.ravel()

        if not gridsize:
            gridsize = length/2

        fig = plt.figure(figsize=(11.7, 8.3))
        plt.subplot(111)
        plt.hexbin(x, y, C=z, gridsize=gridsize, 
            cmap=CM.jet, bins=None)
        plt.axis([x.min(), x.max(), y.min(), y.max()])

        cb = plt.colorbar()
        cb.set_label('Distance')

        return fig

    def get_memberships(self, partition):
        clusters = list(set(partition))
        result = []
        for c in clusters:
            members = []
            for i in range(len(partition)):
                if c == partition[i]:
                    members.append(i)
            result.append(set(members))
        result = sorted(result, key=len, reverse=True)
        return result

    def concatenate_records(self, compound_key, records_list):
        """
        NB: had a version of this which used a reduce construct to concatenate
        the alignments - reduce(lambda x,y: x+y, records_list) - but this led
        to problems of the original object being modified. Deepcopying the 
        first record, ensuring a new memory address for the concatenation, seems 
        more robust.
        """

        partition = self.partitions[compound_key]  # partition = list like [ 1, 1, 2, 1, 1, 2, 3, 3]
        memberships = self.get_memberships(partition)
        clusters = {}
        shorten = lambda k: ''.join([str(k[x])[:3] for x in
                                    range(len(k))])
        index = 1  # index for naming clusters
        for cluster in memberships:
            members = [records_list[i] for i in cluster]
            first = deepcopy(members[0])  # use of deepcopy here is important
            for rec in members[1:]:
                first += rec
            first.name = shorten(compound_key) + '_{0}'.format(index)
            clusters[index] = {'concatenation': first,
                               'members': members}
            index += 1
        self.clusters[compound_key] = Result(clusters)
        return clusters

    def variation_of_information(self, partition_1, partition_2):
        """ 
        Functions to calculate Variation of Information Metric between two 
        clusterings of the same data - SEE Meila, M. (2007). Comparing 
        clusterings: an information based distance. Journal of Multivariate
        Analysis, 98(5), 873-895. doi:10.1016/j.jmva.2006.11.013 

        dependencies:
        math.log
        
        parameters:
        partition_1 (list / array) - a partitioning of a dataset according to 
                some clustering method. Cluster labels are arbitrary.
        partition_2 (list / array) - another partitioning of the same dataset.
                Labels don't need to match, nor do the number of clusters.

        subfunctions:
        get_memberships - parameter partition (list / array)
            returns a list of length equal to the number of clusters found in
            the partition. Each element is the set of members of the cluster.
            Ordering is arbitrary.

        variables used:
            t = total number of points in the dataset
            m1 = cluster memberships from partition_1
            m2 = cluster memberships from partition_2
            l1 = length (i.e. number of clusters) of m1
            l2 = length of m2
            entropy_1 = Shannon entropy of partition_1
            entropy_2 = Shannon entropy of partition_2
            mut_inf = mutual information of partitions
            prob1 = probability distribution of partition 1 - i.e. the 
                probability that a randomly chosen datapoint belongs to
                each cluster (size of cluster / size of dataset)
            prob2 = as above, for partition 2
            intersect = number of common elements in partition 1 [i] and
                partition 2 [j]
        """

        if len(partition_1) != len(partition_2):
            print 'Partition lists are not the same length'
            return 0
        else:
            total = float(len(partition_1))  # Ensure float division later

        m1 = self.get_memberships(partition_1)
        m2 = self.get_memberships(partition_2)
        l1 = len(m1)
        l2 = len(m2)
        entropy_1 = 0
        entropy_2 = 0
        mut_inf = 0
        for i in range(l1):
            prob1 = len(m1[i]) / total
            entropy_1 -= prob1 * log(prob1, 2)
            for j in range(l2):
                if i == 0:  # only calculate these once
                    prob2 = len(m2[j]) / total
                    entropy_2 -= prob2 * log(prob2, 2)
                intersect = len(m1[i] & m2[j])
                if intersect == 0:
                    continue  # because 0 * log(0) = 0 (lim x->0: xlog(x)->0)
                else:
                    mut_inf += intersect / total * log(total
                            * intersect / (len(m1[i]) * len(m2[j])), 2)

        return entropy_1 + entropy_2 - 2 * mut_inf

    # ## Methods for multidimensional scaling

    def get_double_centre(self, matrix):
        """ 
        Double-centres (Gower centres) the input matrix as follows:
        square the input matrix and divide by -2
        from each element subtract the row and column means,
            and add the overall mean
        Returns the double-centred matrix
        """

        matrix = copy(matrix)
        matrix *= matrix
        matrix /= -2.0
        size = len(matrix)
        output = np.zeros([size, size])
        row_means = np.array([np.mean(row) for row in matrix])
        col_means = np.array([np.mean(col) for col in matrix.T])
        col_means.shape = (size, 1)
        matrix_mean = np.mean(matrix)
        matrix -= row_means
        matrix -= col_means
        matrix += matrix_mean
        return matrix

    def get_eigen(self, matrix, standardize=False):
        """
        Calculates the eigenvalues and eigenvectors from the double-
        centred matrix
        Returns a tuple of (eigenvalues, eigenvectors, cumulative
        percentage of variance explained)
        eigenvalues and eigenvectors are sorted in order of eigenvalue
        magnitude, high to low 
        """

        (vals, vecs) = np.linalg.eigh(matrix)
        ind = vals.argsort()[::-1]
        vals = vals[ind]
        vecs = vecs[:, ind]
        cum_var_exp = np.cumsum(100 * abs(vals) / sum(abs(vals)))
        if standardize:
            vecs = vecs * np.sqrt(abs(vals))
        return (vals, vecs, cum_var_exp)

    def get_coords_by_cutoff(
        self,
        vals,
        vecs,
        cum_var_exp,
        cutoff=95,
        normalise=True,
        ):
        """
        Returns fitted coordinates in as many dimensions as are
        needed to explain a given amount of variance (specified 
        in the cutoff)
        """

        i = np.where(cum_var_exp >= cutoff)[0][0]
        coords_matrix = vecs[:, :i + 1]

        if normalise:
            coords_matrix = self.normalise_coords(coords_matrix)
        return coords_matrix

    def get_coords_by_dimension(
        self,
        vals,
        vecs,
        cum_var_exp,
        dimensions=3,
        normalise=True,
        ):
        """
        Returns fitted coordinates in specified number of dimensions,
        and the amount of variance explained)
        """

        coords_matrix = vecs[:, :dimensions]
        varexp = cum_var_exp[dimensions - 1]
        if normalise:
            coords_matrix = self.normalise_coords(coords_matrix)
        return (coords_matrix, varexp)

    def normalise_coords(self, coords_matrix):
        sqsum = np.sum(coords_matrix ** 2,
                       axis=1).reshape(coords_matrix.shape[0], -1)
        return coords_matrix / np.sqrt(sqsum)

    # def affinity_propagation(
    #     self,
    #     distance_matrix,
    #     metric,
    #     n_clusters,
    #     p=None,
    #     ):
    #     """
    #     Uses scikit-learn's AffinityPropagation class to do
    #     affinity propagation clustering.
    #     """

    #     if not metric in self._affinity_propagation:
    #         print '{0} metric not supported by affinity propagation'

    #         # If this message occurs, add missing metric to
    #         # self._affinity_propagation dictionary in __init__

    #     if n_clusters in self._affinity_propagation[metric]:
    #         p = self._affinity_propagation[metric][n_clusters][0]

    #     affinity_matrix = -distance_matrix ** 2
    #     ind = np.triu_indices(len(affinity_matrix), 1)
    #     minval = np.min(affinity_matrix[ind])
    #     medval = np.median(affinity_matrix[ind])
    #     maxval = np.max(affinity_matrix[ind])
    #     step = (medval - minval) / 100

    #     aff = AffinityPropagation(damping=0.75)

    #     if not p:
    #         for pref in np.arange(minval, medval, step):
    #             aff.fit(affinity_matrix, pref)
    #             clusters_found = len(aff.cluster_centers_indices_)
    #             self._affinity_propagation[metric][clusters_found].append(pref)
    #             if clusters_found == n_clusters:
    #                 return aff.labels_

    #         for x in sorted(self._affinity_propagation[metric]):
    #             if x < n_clusters:
    #                 minval = self._affinity_propagation[metric][x][0]
    #             else:
    #                 maxval = self._affinity_propagation[metric][x][0]
    #                 break

    #         step = (maxval - minval) / 200

    #         for pref in np.arange(minval, maxval, step):
    #             aff.fit(affinity_matrix, pref)
    #             clusters_found = len(aff.cluster_centers_indices_)
    #             self._affinity_propagation[metric][clusters_found].append(pref)
    #             if clusters_found == n_clusters:
    #                 return aff.labels_

    #         print 'Unable to find {0} clusters'.format(n_clusters)
    #         return
    #     else:
    #         print 'Using dictionary for {0} clusters: P = {1}'.format(n_clusters,
    #                 p)
    #         aff.fit(affinity_matrix, p)
    #         clusters_found = len(aff.cluster_centers_indices_)
    #         print 'Clusters found = ', clusters_found
    #         assert clusters_found == n_clusters
    #         return aff.labels_

    def spectral(
        self,
        distance_matrix,
        prune=True,
        sigma7=False,
        ):
        """
        1st: Calculates an affinity matrix from a distance matrix, using the
        local scaling transform from Zelnik-Manor and Perona (2004):
        affinity[i,j] = exp(-distance[i,j]^2/(sigma_i*sigma_j)).
        2nd: Returns a normalised Laplacian matrix from the affinity matrix.
        Optionally the similarity matrix can be pruned using a k-nearest neighbours
        approach as in Leigh et al. (2011).
        Note: the normalised Laplacian according to Ng et al. is different to the
        normalised Laplacian according to Von Luxburg: 
        Ng et al. give D(-1/2). W. D(-1/2)
        VL gives I - D(-1/2). W. D(-1/2)

        References: 
        Luxburg, U. (2007). A tutorial on spectral clustering. 
        Statistics and Computing, 17(4), 395-416. 
        doi:10.1007/s11222-007-9033-z

        Leigh, J. W. et al. (2011). 
        Let Them Fall Where They May: Congruence Analysis in Massive 
        Phylogenetically Messy Data Sets. 
        Molecular Biology and Evolution, 28(10), 2773-2785.
        doi:10.1093/molbev/msr110

        P Perona and L. Zelnik-Manor. (2004).
        Self-tuning spectral clustering.
        Advances in neural information processing systems, 2004 vol. 17 pp. 1601-1608
        """

        size = len(distance_matrix)  # assume square and symmetrical input

        def isconnected(matrix):
            """
            Checks that all nodes are reachable from the first node - i.e. that
            the graph is fully connected. The approach is borrowed from 
            isconnected function from graph.c in Leigh's Conclustador program.
            """

            # INIT

            matrix = np.array(matrix)
            num_nodes = len(matrix)
            checklist = [0]
            checklength = 1
            reachable = [0]

            # ALGORITHM

            while checklength > 0:
                node = checklist.pop(0)
                checklength -= 1
                for edge in range(num_nodes):
                    if matrix[node, edge] == 1:
                        reachable_node = edge
                        if reachable_node not in reachable:
                            reachable.append(reachable_node)
                            checklist.append(reachable_node)
                            checklength += 1

            # RESULT

            for i in range(num_nodes):
                if not i in reachable:
                    return False
            return True

        def nodivzero(d):
            if 0 in d.values():
                return False
            else:
                return True

        def knn(matrix, k):
            """
            Acts on distance matrix. For each datapoint, finds
            the `k` nearest neighbours. Returns an adjacency
            matrix, and a dictionary of the kth distance for 
            each node.
            """

            kneighbour_matrix = np.zeros([size, size])
            max_dists = {}
            for i in range(size):
                sorted_dists = matrix[i].argsort()
                for j in sorted_dists[:k]:
                    kneighbour_matrix[i, j] = kneighbour_matrix[j, i] = \
                        1
                    max_dists[i] = matrix[i, sorted_dists[k - 1]]
            return (kneighbour_matrix, max_dists)

        def get_affinity_matrix(distance_matrix, kneighbour_matrix,
                                max_dists):
            """
            Makes weighted adjacency matrix along the lines of
            Zelnik-Manor and Perona (2004), with local scaling.
            """

            affinity_matrix = np.zeros([size, size])
            for i in range(size):
                for j in range(size):
                    if i != j and kneighbour_matrix[i, j] == 1:
                        distance = distance_matrix[i, j]
                        sigma_i = max_dists[i]
                        sigma_j = max_dists[j]
                        affinity_matrix[i, j] = np.exp(-distance ** 2
                                / (sigma_i * sigma_j))
            return affinity_matrix

        def laplace(affinity_matrix):

            D = np.diag(affinity_matrix.sum(axis=1))
            invRootD = np.sqrt(np.linalg.inv(D))
            return np.dot(np.dot(invRootD, affinity_matrix), invRootD)

        if prune:  # 'guess a number' strategy
            mink = 1
            maxk = size
            guessk = int(np.log(size).round())
            while maxk - mink != 1:
                test = knn(distance_matrix, guessk)
                if isconnected(test[0]) and nodivzero(test[1]):

                    # either correct or too high
                    # try a lower number

                    maxk = guessk
                    guessk = mink + (guessk - mink) / 2
                else:

                    # too low

                    mink = guessk
                    guessk = guessk + (maxk - guessk) / 2
            (kneighbour_matrix, max_dists) = knn(distance_matrix,
                    guessk + 1)
        else:

            (kneighbour_matrix, max_dists) = knn(distance_matrix, size)

        affinity_matrix = get_affinity_matrix(distance_matrix,
                kneighbour_matrix, max_dists)
        md7 = knn(distance_matrix, 7)[1]
        if nodivzero(md7):
            am7 = get_affinity_matrix(distance_matrix,
                    kneighbour_matrix, md7)
        else:
            print 'Setting sigma(i) > d(S(i),S(7)) to avoid dividing by zero.'
            if prune:
                print 'Sigma = d(S(i),S({0})'.format(guessk + 1)
            else:
                print 'Sigma = d(S(i),S({0})'.format(size)
            sigma7 = False

        if sigma7:
            return laplace(am7)
        else:
            return laplace(affinity_matrix)
