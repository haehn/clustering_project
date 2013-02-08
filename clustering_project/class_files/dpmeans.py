#!/usr/bin/env python


class DP_tree(object):

    (sumD_old, sumD_new) = (np.Inf, np.Inf)

    def __init__(self, data=None, lmbda=None):
        self.data = data
        self.lmbda = lmbda

    def __str__(self):
        pass

    def add_data(self, data):
        self.data = data
        self.assignments = np.array([0 for _ in data])

    def distance(self, tree1, tree2):
        pass

    def cluster_centre(self, cluster):
        pass

    def members(self):
        pass

    def converged(self):
        pass
