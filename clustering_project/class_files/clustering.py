#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import dendropy as dpy
import GeoMeTreeHack
from result import Result
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from Bio.Cluster import kmedoids
import matplotlib.pyplot as plt
from math import log
import copy


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
        if metric == 'geodesic':
            geotrees = [tree.newick for tree in trees]
            for i in range(num_trees):
                for j in range(i + 1, num_trees):
                    matrix[i][j] = matrix[j][i] = \
                        GeoMeTreeHack.main(geotrees[i], geotrees[j])
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
                        matrix[i][j] = matrix[j][i] = \
                            dpytrees[i].symmetric_difference(dpytrees[j])
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
        ):
        matrix = self.get_distance_matrix(trees, metric, invert=invert,\
            normalise=normalise)
        self.distance_matrices[metric] = matrix

    def order(self, l, num=1):
        """ The clustering returned by the hcluster module gives 
            group membership without regard for numerical order 
            This function preserves the group membership, but sorts 
            the labelling into numerical order
        """

        # base case

        if num >= max(l):
            return l
        else:

        # recursion on num

            outl = []
            change_places = None
            for i in range(len(l)):
                if l[i] < num:
                    outl.append(l[i])
                else:
                    change_places = l[i]
                    break
            for j in range(i, len(l)):
                if l[j] == change_places:
                    outl.append(num)
                elif l[j] == num:
                    outl.append(change_places)
                else:
                    outl.append(l[j])
            return self.order(outl, num + 1)

    def put_partition(
        self,
        metric,
        linkage_method,
        nclasses,
        names,
        criterion='distance',
        ):
        """ Returns list of cluster assignments from linkage
            Criteria: 'maxclust' -  set threshold as (maximum) number 
                                    of groups to cluster into 
                      'distance' -  set threshold as cutpoint at which 
                                    to separate dendrogram into clusters 
                                    (threshold in range float(0,1)) """

        dm = self.distance_matrices[metric]
        compound_key = (metric, linkage_method, nclasses)

        if linkage_method == 'kmedoids':
            p = []
            for i in range(100):
                p.append(kmedoids(dm, nclusters=nclasses, npass=100))
            T = sorted(p, key=lambda x:x[1])[0][0]

        else:
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
            self.plotting_info[compound_key] = (linkmat, names, threshold)
        
        T = self.order(T)  # puts cluster labels in ascending order
        self.partitions[compound_key] = T
        
    def plot_dendrogram(self, compound_key):
        """
        Extracts data from clustering to plot dendrogram
        """

        partition = self.partitions[compound_key]
        (linkmat, names, threshold) = self.plotting_info[compound_key]
        fig = plt.figure(figsize=(11.7,8.3))
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
        partition = self.partitions[compound_key] # partition = list like [ 1, 1, 2, 1, 1, 2, 3, 3]
        memberships = self.get_memberships(partition)
        clusters = {}
        shorten = lambda k: ''.join([str(k[x])[:3] for x in
                                    range(len(k))])
        index = 1  # index for naming clusters
        for cluster in memberships:
            members = [records_list[i] for i in cluster]
            first = copy.deepcopy(members[0]) # use of deepcopy here is important
            for rec in members[1:]:
                first += rec
            first.name = shorten(compound_key) + '_{0}'.format(index)
            clusters[index] = {'concatenation': first, 'members': members}
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
