#!/usr/bin/env python

import_debugging = False
if import_debugging:
    print 'clustering.py imports:'
import numpy as np
if import_debugging:
    print '  numpy (cl)'
np.set_printoptions(threshold='nan', precision=3)
import sys
if import_debugging:
    print '  sys (cl)'
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
if import_debugging:
    print '  scipy.cluster.hierarchy::linkage, fcluster, dendrogram (cl)'
from Bio.Cluster import kmedoids
if import_debugging:
    print '  Bio.Cluster::kmedoids (cl)'
import matplotlib.pyplot as plt
if import_debugging:
    print '  matplotlib.pyplot (cl)'
from sklearn.cluster import KMeans
if import_debugging:
    print '  sklearn.cluster::KMeans (cl)'
from collections import defaultdict
if import_debugging:
    print '  collections::defaultdict (cl)'
import evrot
if import_debugging:
    print '  evrot (cl)'


class Clustering(object):

    """
    Apply clustering methods to distance matrix
    """

    def __init__(self):
        """
        self.partitions uses a compound key to retrieve partitions
        key = tuple of (distance_metric, linkage_method, num_classes)
        """

        self.cache = {}

    def __str__(self):
        pass

        # s = ''
        # num_partitions = len(self.partitions)
        # num_clusters = len(self.clusters)
        # if num_clusters == 1:
        #     s += '{0} partition calculated:\n'.format(num_partitions)
        # else:
        #     s += '{0} partitions calculated:\n'.format(num_partitions)
        # for p in self.partitions:
        #     s += ' '.join(str(x) for x in p) + '\n'
        # return s

    def clear_cache(self):
        self.cache = {}

    def run_kmedoids(self, dm, nclusters):

        if dm.metric == 'rf':
            matrix = dm.add_noise(dm.matrix)
        else:
            matrix = dm.matrix

        p = [kmedoids(matrix, nclusters=nclusters, npass=100) for _ in
             range(100)]
        p.sort(key=lambda x: x[1])
        T = self.order(p[0][0])
        return T

    def run_spectral(
        self,
        dm,
        nclusters,
        prune=True,
        sigma7=False,
        recalculate=False,
        ):

        if dm.metric == 'rf':
            noise = True
        else:
            noise = False

        if recalculate or not 'spectral_decomp' in self.cache:
            laplacian = self.spectral(dm, prune=prune, add_noise=noise)

            (eigvals, eigvecs, cve) = self.get_eigen(laplacian,
                    standardize=False)
            self.cache['spectral_decomp'] = (eigvals, eigvecs, cve)
            self.cache['laplacian'] = laplacian
        else:

            (eigvals, eigvecs, cve) = self.cache['spectral_decomp']

        coords = self.get_coords_by_dimension(eigvals, eigvecs, cve,
                nclusters, normalise=True)[0]
        T = self.run_KMeans(coords, nclusters)
        return T

    def run_spectral_rotate(
        self,
        dm,
        prune=True,
        KMeans=True,
        recalculate=False,
        max_groups=None,
        min_groups=2,
        verbose=True,
        ):

        if dm.metric == 'rf':
            noise = True
        else:
            noise = False
        if recalculate or not 'spectral_decomp' in self.cache:
            laplacian = self.spectral(dm, prune=prune, add_noise=noise)

            (eigvals, eigvecs, cve) = self.get_eigen(laplacian,
                    standardize=False)
            self.cache['spectral_decomp'] = (eigvals, eigvecs, cve)
            self.cache['laplacian'] = laplacian
        else:

            (eigvals, eigvecs, cve) = self.cache['spectral_decomp']

        # ######################
        # CLUSTER_ROTATE STUFF HERE

        M = dm.matrix
        if not max_groups:
            max_groups = int(np.sqrt(M.shape[0]) + np.power(M.shape[0],
                             1.0 / 3))
        (nclusters, clustering, quality_scores, rotated_vectors) = \
            self.cluster_rotate(eigvecs, max_groups=max_groups,
                                min_groups=min_groups)

        translate_clustering = [None] * len(dm.M)
        no_of_empty_clusters = 0
        for (group_number, group_membership) in enumerate(clustering):
            if len(group_membership) == 0:
                no_of_empty_clusters += 1
            for index in group_membership:
                translate_clustering[index - 1] = group_number
        clustering = self.order(translate_clustering)
        if no_of_empty_clusters > 0:
            print 'Subtracting {0} empty {1}'.format(no_of_empty_clusters,
                    ('cluster' if no_of_empty_clusters
                    == 1 else 'clusters'))
            nclusters -= no_of_empty_clusters

        # ######################

        if verbose:
            print 'Discovered {0} clusters'.format(nclusters)
            print 'Quality scores: {0}'.format(quality_scores)
            if KMeans:
                print 'Pre-KMeans clustering: {0}'.format(clustering)
        if KMeans:
            T = self.run_KMeans(rotated_vectors, nclusters=nclusters)
        else:
            T = clustering
        return (T, nclusters, quality_scores)

        # ######################

    def run_NJW(
        self,
        dm,
        nclusters,
        recalculate=False,
        ):

        if dm.metric == 'rf':
            noise = True
        else:
            noise = False
        if recalculate or not 'NJW_decomp' in self.cache:
            laplacian = self.NJW(matrix, sigma=np.median(matrix),
                                 add_noise=noise)

            (eigvals, eigvecs, cve) = self.get_eigen(laplacian,
                    standardize=False)
            self.cache['NJW_decomp'] = (eigvals, eigvecs, cve)
            self.cache['NJW_laplacian'] = laplacian
        else:

            (eigvals, eigvecs, cve) = self.cache['NJW_decomp']

        coords = self.get_coords_by_dimension(eigvals, eigvecs, cve,
                nclusters, normalise=True)[0]
        T = self.run_KMeans(coords, nclusters)
        return T

    def run_ShiMalik(
        self,
        dm,
        nclusters,
        recalculate=False,
        ):

        if dm.metric == 'rf':
            noise = True
        else:
            noise = False
        if recalculate or not 'SM_decomp' in self.cache:
            laplacian = self.ShiMalik(matrix, sigma=np.median(matrix),
                    add_noise=noise)

            (eigvals, eigvecs, cve) = self.get_eigen(laplacian,
                    standardize=False)
            self.cache['SM_decomp'] = (eigvals, eigvecs, cve)
            self.cache['SM_laplacian'] = laplacian
        else:

            (eigvals, eigvecs, cve) = self.cache['SM_decomp']

        coords = self.get_coords_by_dimension(eigvals, eigvecs, cve,
                nclusters, normalise=True)[0]
        T = self.run_KMeans(coords, nclusters)
        return T

    def run_hierarchical(
        self,
        dm,
        nclusters,
        linkage_method,
        ):

        if dm.metric == 'rf':
            matrix = dm.add_noise(dm.matrix)
        else:
            matrix = dm.matrix

        linkmat = linkage(matrix, linkage_method)
        linkmat_size = len(linkmat)
        if nclusters <= 1:
            br_top = linkmat[linkmat_size - nclusters][2]
        else:
            br_top = linkmat[linkmat_size - nclusters + 1][2]
        if nclusters >= len(linkmat):
            br_bottom = 0
        else:
            br_bottom = linkmat[linkmat_size - nclusters][2]
        threshold = 0.5 * (br_top + br_bottom)
        T = fcluster(linkmat, threshold, criterion='distance')
        T = self.order(T)
        return T

    def run_MDS(
        self,
        dm,
        nclusters,
        recalculate=False,
        ):

        if recalculate or not 'MDS_decomp' in self.cache:

            if dm.metric == 'rf':
                noise = True
            else:
                noise = False

            dbc = dm.get_double_centre(add_noise=noise)
            (eigvals, eigvecs, cve) = self.get_eigen(dbc,
                    standardize=True)
            self.cache['MDS_decomp'] = (eigvals, eigvecs, cve)
            self.cache['dbc'] = dbc
        else:

            (eigvals, eigvecs, cve) = self.cache['MDS_decomp']

        coords = self.get_coords_by_cutoff(eigvals, eigvecs, cve, 95,
                normalise=False)
        T = self.run_KMeans(coords, nclusters)
        return T

    def run_KMeans(self, coords, nclusters):
        est = KMeans(n_clusters=nclusters)
        est.fit(coords)
        T = self.order(est.labels_)
        return T

    def run_clustering(
        self,
        dm,
        method,
        nclusters,
        prune=True,
        recalculate=False,
        ):

        if method == 'kmedoids':
            return self.run_kmedoids(dm, nclusters)
        elif method == 'spectral':
            return self.run_spectral(dm, nclusters, prune,
                    recalculate=recalculate)
        elif method == 'spectral_rotate':
            return self.run_spectral_rotate(dm, recalculate=recalculate)
        elif method == 'spectral-prune':
            return self.run_spectral(dm, nclusters, prune=False,
                    sigma7=True, recalculate=recalculate)
        elif method == 'NJW':
            return self.run_NJW(dm, nclusters, recalculate=recalculate)
        elif method == 'ShiMalik':
            return self.run_ShiMalik(dm, nclusters,
                    recalculate=recalculate)
        elif method == 'MDS':
            return self.run_MDS(dm, nclusters, recalculate=recalculate)
        elif method in ['single', 'complete', 'average', 'ward']:
            return self.run_hierarchical(dm, nclusters, method)
        else:
            print 'Unrecognised method: {0}'.format(method)

    def order(self, l):
        """
        The clustering returned by the hcluster module gives 
        group membership without regard for numerical order 
        This function preserves the group membership, but sorts 
        the labelling into numerical order
        """

        list_length = len(l)

        d = defaultdict(list)
        for (i, element) in enumerate(l):
            d[element].append(i)

        l2 = [None] * list_length

        for (name, index_list) in enumerate(sorted(d.values(),
                key=min), start=1):
            for index in index_list:
                l2[index] = name

        return tuple(l2)

    def plot_embedding(
        self,
        dm,
        metric,
        linkage,
        nclasses,
        dimensions=3,
        embed='MDS',
        standardize=False,
        normalise=True,
        ):

        if not dimensions in [2, 3]:
            print '2D or 3D only'
            return

        dm = self.distance_matrices[metric]
        partition = self.partitions[(metric, linkage, nclasses)]

        if embed == 'MDS':
            dbc = dm.get_double_centre()
            (eigvals, eigvecs, cve) = self.get_eigen(dbc,
                    standardize=standardize)
            (coords, varexp) = self.get_coords_by_dimension(eigvals,
                    eigvecs, cve, dimensions, normalise=normalise)
        elif embed == 'spectral':

            laplacian = self.spectral(dm)
            (eigvals, eigvecs, cve) = self.get_eigen(laplacian,
                    standardize=standardize)
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

    # ## Methods for eigen decomposition

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

    def spectral(
        self,
        dm,
        prune=True,
        sigma7=False,
        add_noise=False,
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

        size = len(dm.matrix)  # assume square and symmetrical input
        if size <= 10:  # no point pruning a small matrix
            prune = False

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

        # def knn(dm, k, add_noise):
        #     """
        #     Acts on distance matrix. For each datapoint, finds
        #     the `k` nearest neighbours. Returns an adjacency
        #     matrix, and a dictionary of the kth distance for 
        #     each node.
        #     """

        #     return dm.get_knn(k)

        # def get_affinity_matrix(dm, kneighbour_matrix, max_dists):
        #     """
        #     Makes weighted adjacency matrix along the lines of
        #     Zelnik-Manor and Perona (2004), with local scaling.
        #     """

        #     return dm.get_affinity_matrix(dm, kneighbour_matrix,
        #             max_dists)

        def laplace(affinity_matrix):

            diagonal = affinity_matrix.sum(axis=1) \
                - affinity_matrix.diagonal()
            invRootD = np.diag(np.sqrt(1 / diagonal))
            return np.dot(np.dot(invRootD, affinity_matrix), invRootD)

        # prune the adjacency matrix

        if prune:  # binary search
            mink = 1
            maxk = size
            guessk = int(np.log(size).round())
            while maxk - mink != 1:
                test = dm.get_knn(guessk, add_noise=add_noise)
                if isconnected(test[0]) and nodivzero(test[1]):

                    # either correct or too high
                    # try a lower number

                    maxk = guessk
                    guessk = mink + (guessk - mink) / 2
                else:

                    # too low

                    mink = guessk
                    guessk = guessk + (maxk - guessk) / 2
            (kneighbour_matrix, max_dists) = \
                dm.get_knn(guessk + 1,
                           add_noise=add_noise)
            print 'Pruning adjacencies to {0}-NN'.format(guessk + 1)
        else:

            (kneighbour_matrix, max_dists) = \
                dm.get_knn(size, add_noise=add_noise)

        affinity_matrix = dm.get_affinity_matrix(kneighbour_matrix,
                max_dists, add_noise=add_noise)

        # Tune the sigma parameter

        md7 = dm.get_knn(7, add_noise=False)[1]
        if nodivzero(md7):
            affinity_matrix = dm.get_affinity_matrix(kneighbour_matrix,
                    md7, add_noise=add_noise)
            print 'Setting sigma to 7-NN'
        else:

            if prune and nodivzero(max_dists):
                print 'Setting sigma to {0}-NN'.format(guessk + 1)
            else:
                print 'Setting sigma to {0}-NN'.format(size)
                mdmax = dm.get_knn(size,
                                   add_noise=add_noise)[1]
                affinity_matrix = \
                    dm.get_affinity_matrix(kneighbour_matrix, mdmax,
                        add_noise=add_noise)

        # normalise affinities to [0,1]

        affinity_matrix *= 1.0 / np.max(affinity_matrix)
        L = laplace(affinity_matrix)
        return L

    def NJW(self, distance_matrix, sigma):
        size = distance_matrix.shape[0]
        A = np.exp(-distance_matrix ** 2 / (2 * sigma))
        A.flat[::size + 1] = 0.0
        D = np.diag(A.sum(axis=1))
        invRootD = np.sqrt(np.linalg.inv(D))
        L = np.dot(invRootD, np.dot(A, invRootD))
        return L

    def ShiMalik(self, distance_matrix, sigma):
        size = distance_matrix.shape[0]
        A = np.exp(-distance_matrix ** 2 / (2 * sigma))
        A.flat[::size + 1] = 0.0
        D = np.diag(A.sum(axis=1))
        invD = np.linalg.inv(D)
        L = np.dot(invD, A)
        return L

    def cluster_rotate(
        self,
        eigenvectors,
        max_groups,
        min_groups=2,
        ):

        groups = range(min_groups, max_groups + 1)
        vector_length = eigenvectors.shape[0]
        current_vector = eigenvectors[:, :groups[0]]
        n = max_groups - min_groups + 1

        quality_scores = [None] * n
        clusters = [None] * n
        rotated_vectors = [None] * n

        for g in range(n):
            if g > 0:
                current_vector = np.concatenate((rotated_vectors[g
                        - 1], eigenvectors[:, groups[g] - 1:
                        groups[g]]), axis=1)

            (clusters[g], quality_scores[g], rotated_vectors[g]) = \
                evrot.main(current_vector)

        # Find the highest index of quality scores where the
        # score is within 0.0025 of the maximum:
        # this is our chosen number of groups

        max_score = max(quality_scores)
        index = quality_scores.index(max_score)
        start = index + 1
        for (i, score) in enumerate(quality_scores[index + 1:],
                                    start=start):
            if abs(score - max_score) < 0.0025:
                index = i

        return (groups[index], clusters[index], quality_scores,
                rotated_vectors[index])
