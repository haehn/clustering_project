#!/usr/bin/env python

import numpy as np
import os
import dendropy as dpy
from matplotlib import pyplot as plt
from matplotlib import cm as CM


class DistanceMatrix(object):

    metrics_dict = {
        'rf': 'Robinson-Foulds distances',
        'wrf': 'Weighted Robinson-Foulds distances',
        'euc': 'Euclidean distances',
        'geo': 'Geodesic distances',
        None: 'Empty matrix',
        }

    def __init__(
        self,
        trees,
        tmpdir='/tmp',
        gtp_path='./class_files',
        ):

        size = len(trees)
        self.matrix = np.zeros((size, size), dtype='float')
        self.metric = None
        self.trees = trees
        self.tmpdir = tmpdir
        self.gtp_path = gtp_path

    def __str__(self):
        return '\n'.join([str(self.matrix),
                         'Metric: {0}'.format(self.metrics_dict[self.metric])])

    def convert_to_dendropy_trees(self, trees=None):
        taxa = dpy.TaxonSet()
        if not trees:
            trees = self.trees
        dpy_tree_list = [dpy.Tree.get_from_string(tree.newick, 'newick'
                         , taxon_set=taxa) for tree in trees]
        return dpy_tree_list

    def get_rf_distance(self, tree1, tree2):
        return tree1.symmetric_difference(tree2)

    def get_wrf_distance(self, tree1, tree2):
        return tree1.robinson_foulds_distance(tree2)

    def get_euc_distance(self, tree1, tree2):
        return tree1.euclidean_distance(tree2)

    def get_distance_matrix(
        self,
        metric,
        trees=None,
        normalise=False,
        tmpdir=None,
        gtp_path=None,
        ):
        """
        Generates pairwise distance matrix between trees
        Uses one of the following distance metrics:
        Robinson-Foulds distance - topology only (='rf')
        Robinson-Foulds distance - branch lengths (='wrf')
        Euclidean distance - Felsenstein's branch lengths
            distance (='euc')
        Geodesic distance - branch lengths (='geo')
        """

        matrix = self.matrix
        if not trees:
            trees = self.trees
        if not tmpdir:
            tmpdir = self.tmpdir
        if not gtp_path:
            gtp_path = self.gtp_path
        assert os.path.isfile('{0}/gtp.jar'.format(gtp_path))
        num_trees = matrix.shape[0]
        self.metric = metric

        if metric == 'geo':
            rooted = all([tree.rooted for tree in trees])
            with open('{0}/geotrees.nwk'.format(tmpdir), 'w') as file:
                file.write('\n'.join([tree.newick.rstrip() for tree in
                           trees]))
            if rooted:
                os.system('java -jar {1}/gtp.jar -o {0}/output.txt {0}/geotrees.nwk'.format(tmpdir,
                          gtp_path))
            else:
                os.system('java -jar {1}/gtp.jar -u -o {0}/output.txt {0}/geotrees.nwk'.format(tmpdir,
                          gtp_path))
            with open('{0}/output.txt'.format(tmpdir)) as file:
                for line in file:
                    line = line.rstrip()
                    if line:
                        (i, j, value) = line.split()
                        i = int(i)
                        j = int(j)
                        value = float(value)
                        matrix[i, j] = matrix[j, i] = value
            os.remove('{0}/output.txt'.format(tmpdir))
            os.remove('{0}/geotrees.nwk'.format(tmpdir))
        elif metric == 'rf':

            dpy_trees = self.convert_to_dendropy_trees(trees)
            ntax = len(dpy_trees[0].leaf_nodes())
            max_rf = 2.0*(ntax-3)
            for i in range(num_trees):
                for j in range(i + 1, num_trees):
                    distance = self.get_rf_distance(dpy_trees[i],
                            dpy_trees[j])
                    distance /= max_rf
                    matrix[i, j] = matrix[j, i] = distance
        elif metric == 'wrf':

            dpy_trees = self.convert_to_dendropy_trees(trees)
            for i in range(num_trees):
                for j in range(i + 1, num_trees):
                    distance = self.get_wrf_distance(dpy_trees[i],
                            dpy_trees[j])
                    matrix[i, j] = matrix[j, i] = distance
        elif metric == 'euc':

            dpy_trees = self.convert_to_dendropy_trees(trees)
            for i in range(num_trees):
                for j in range(i + 1, num_trees):
                    distance = self.get_euc_distance(dpy_trees[i],
                            dpy_trees[j])
                    matrix[i, j] = matrix[j, i] = distance
        else:

            print 'Unrecognised distance metric'
            return

        if normalise:
            matrix = matrix / np.max(matrix)

        return matrix

    def add_noise(self, dm=None):
        if dm is None:
            dm = self.matrix
        size = len(dm)
        new = np.zeros((size, size))
        r = range(size)
        for i in r:
            for j in r[i+1:]:
                eps = np.random.normal(0, 0.001)
                if dm[i, j] + eps > 0:
                    new[i, j] = new[j, i] = dm[i, j] + eps
                else:
                    new[i, j] = new[j, i] = dm[i, j] - eps
        return new

    def plot_heatmap(self, sort_partition=None):
        """
        Sort partition should be a flatlist of the clusters
        as returned by Partition().get_memberships(..., flatten=True)
        """

        dm = np.array(self.matrix, copy=True)
        length = dm.shape[0]
        datamax = abs(dm).max()
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

    def check_euclidean(self):
        """
        A distance matrix is euclidean iff
        F = -0.5 * (I - 1/n)D(I - 1/n) is PSD,
        where I is the identity matrix
              D is the distance matrix
              1 is a square matrix of ones
              n is the matrix size, common to all
        """

        D = self.matrix
        I = np.identity( D.shape )
        ones = np.ones( D.shape )

        bracket = I - ones/D.shape[0]
        F = -0.5 * bracket.dot(D.dot(bracket))

        # Can calculate Cholesky decomp iff F is PSD
        try:
            np.linalg.cholesky(F)
        except: return False
        
        return True