#!/usr/bin/env python

from matplotlib import pyplot as plt
from matplotlib import cm as CM
import numpy as np
from externals import GTP
from utils.dpy import convert_to_dendropy_trees, get_euc_distance, \
    get_rf_distance, get_wrf_distance


class DistanceMatrix(object):

    def __init__(
        self,
        trees,
        metric,
        tmpdir='/tmp',
        **kwargs
        ):
        """ Initialise with a list of trees to compare, and a metric (one of
        'geo', 'euc', 'wrf', 'rf') to compare them """

        self.tmpdir = tmpdir
        self.matrix = self.get_distance_matrix(trees, metric)
        self.metric = metric

    def __str__(self):
        return '\n'.join([str(self.matrix),
                         'Metric: {0}'.format(self.metrics_dict[self.metric])])

    def copy(self, add_noise=False):
        copy = self.__new__(type(self))
        copy.__dict__ = dict(self.__dict__.items())
        if add_noise:
            return MatrixTricks(copy).add_noise()
        return copy

    def get_dendropy_distances(self, dpy_trees, fn):
        num_trees = len(dpy_trees)

        matrix = np.zeros((num_trees, num_trees))
        for i in range(num_trees):
            for j in range(i + 1, num_trees):
                distance = fn(dpy_trees[i], dpy_trees[j])
                matrix[i, j] = matrix[j, i] = distance
        return matrix

    def get_geo_distances(self, trees, tmpdir=None):

        tmpdir = tmpdir or self.tmpdir

        g = GTP(tmpdir=tmpdir)
        return g.run(trees)

    def get_distance_matrix(
        self,
        trees,
        metric,
        tmpdir=None,
        ):
        """ Generates pairwise distance matrix between trees Uses one of the
        following distance metrics: Robinson-Foulds distance - topology only
        (='rf') Robinson-Foulds distance - branch lengths (='wrf') Euclidean
        distance - Felsenstein's branch lengths distance (='euc') Geodesic
        distance - branch lengths (='geo') """

        tmpdir = tmpdir or self.tmpdir
        if metric == 'geo':
            return self.get_geo_distances(trees, tmpdir=tmpdir)

        dpy_trees = convert_to_dendropy_trees(trees)
        
        if metric == 'rf':
            matrix = self.get_dendropy_distances(dpy_trees, get_rf_distance)
        
        elif metric == 'wrf':
            matrix = self.get_dendropy_distances(dpy_trees, get_wrf_distance)
        
        elif metric == 'euc':
            matrix = self.get_dendropy_distances(dpy_trees, get_euc_distance)
        
        else:
            print 'Unrecognised distance metric'
            return

        return matrix


class MatrixTricks(object):

    def __init__(self, dm):

        self.matrix = dm.matrix
        self.shape = dm.matrix.shape
        self.size = dm.matrix.size

    def __str__(self):
        return '\n'.join([str(self.matrix),
                         'Metric: {0}'.format(self.metrics_dict[self.metric])])

    def add_noise(self):

        noise = np.random.normal(0, 0.001, self.size).reshape(self.shape)
        np.fill_diagonal(noise, 0.0)
        noisy = self.matrix + noise
        noisy[noisy < 0] = np.abs(noisy[noisy < 0])
        return noisy

    def affinity_matrix(
        self,
        kneighbour_matrix,
        max_dists,
        local_scaling=True,
        sigma=2,
        add_noise=False,
        ):
        """ Makes weighted adjacency matrix along the lines of Zelnik-Manor and
        Perona (2004), with local scaling. """

        if add_noise:
            M = self.add_noise()
        else:
            M = self.matrix
        shape = M.shape
        affinity_matrix = np.zeros(shape)
        for i in range(shape[0]):
            for j in range(shape[1]):
                if kneighbour_matrix[i, j] == 1:
                    distance = M[i, j]
                    if local_scaling:
                        sigma_i = max_dists[i]
                        sigma_j = max_dists[j]
                        affinity_matrix[i, j] = np.exp(-distance ** 2
                                / (sigma_i * sigma_j))
                    else:
                        affinity_matrix[i, j] = np.exp(-distance ** 2 / 2
                                * sigma)
        return affinity_matrix

    def check_euclidean(self):
        """ A distance matrix is euclidean iff F = -0.5 * (I - 1/n)D(I - 1/n) is
        PSD, where I is the identity matrix D is the distance matrix 1 is a
        square matrix of ones n is the matrix size, common to all """

        D = self.matrix
        I = np.identity(D.shape)
        ones = np.ones(D.shape)

        bracket = I - ones / D.shape[0]
        F = -0.5 * bracket.dot(D.dot(bracket))

        try:
            np.linalg.cholesky(F) # Can calculate Cholesky decomp iff F is PSD
        except:
            return False

        return True

    def double_centre(self, square_input=False, add_noise=False):
        """ Double-centres the input matrix: From each element: Subtract the row
        mean Subtract the column mean Add the grand mean Divide by -2 Method
        from: Torgerson, W S (1952). Multidimensional scaling: I. Theory and
        method. """

        if add_noise:
            M = self.add_noise()
        else:
            M = np.array(self.matrix, copy=True)
        if square_input:
            M *= M
        (rows, cols) = M.shape
        cm = np.mean(M, axis=0)  # column means
        rm = np.mean(M, axis=1).reshape((rows, 1))  # row means
        gm = np.mean(cm)  # grand mean
        M -= rm + cm - gm
        M /= -2
        return M

    def eigen(self, matrix, standardize=False):
        """ Calculates the eigenvalues and eigenvectors from the double-centred
        matrix. Returns a tuple of (eigenvalues, eigenvectors, cumulative
        percentage of variance explained). Eigenvalues and eigenvectors are
        sorted in order of eigenvalue magnitude, high to low """

        (vals, vecs) = np.linalg.eigh(matrix)
        ind = vals.argsort()[::-1]
        vals = vals[ind]
        vecs = vecs[:, ind]
        cum_var_exp = np.cumsum(100 * abs(vals) / sum(abs(vals)))
        if standardize:
            vecs = vecs * np.sqrt(abs(vals))
        return (vals, vecs, cum_var_exp)

    def knn(self, k, add_noise=False):
        """ Acts on distance matrix. For each datapoint, finds the `k` nearest
        neighbours. Returns an adjacency matrix, and a dictionary of the kth
        distance for each node. """

        if add_noise:
            M = self.add_noise()
        else:
            M = self.matrix
        shape = M.shape
        kneighbour_matrix = np.zeros(shape)
        max_dists = {}
        for i in range(shape[0]):
            sorted_dists = M[i].argsort()
            for j in sorted_dists[:k]:
                kneighbour_matrix[i, j] = kneighbour_matrix[j, i] = 1
                max_dists[i] = M[i, sorted_dists[k - 1]]
        return (kneighbour_matrix, max_dists)


class MatrixPlots(object):

    def __init__(self, dm):

        self.matrix = dm.matrix.copy()
        self.shape = dm.matrix.shape
        self.size = dm.matrix.size

    def get_permutation_matrix(self, input_ordering, desired_ordering):
        length = len(input_ordering)
        if not len(desired_ordering) == length:
            print 'List lengths don\'t match'
            return
        P = np.zeros((length, length), dtype=np.int)
        for i in range(length):
            j = desired_ordering.index(input_ordering[i])
            P[i, j] = 1
        return P

    def plot_heatmap(self, sort_partition=None):
        """ Sort partition should be a flatlist of the clusters as returned by
        Partition().get_memberships(..., flatten=True) """

        dm = self.matrix
        length = self.shape[0]
        datamax = np.abs(dm).max()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ticks_at = [0, 0.5 * datamax, datamax]
        if sort_partition:
            p = self.get_permutation_matrix(range(len(sort_partition)),
                    sort_partition)
            dm = np.dot(p.T, np.dot(dm, p))
        cax = ax.imshow(
            dm,
            interpolation='nearest',
            origin='lower',
            extent=[0.0, length, 0.0, length],
            vmin=0,
            vmax=datamax,
            cmap=CM.Blues,
            )
        cbar = fig.colorbar(cax, ticks=ticks_at, format='%1.2g')
        cbar.set_label('Distance')
        return fig
