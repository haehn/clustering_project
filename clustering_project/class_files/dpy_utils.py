#!/usr/bin/env python

import dendropy as dpy

# DENDROPY UTILS

def get_rf_distance(tree1, tree2):
    return tree1.symmetric_difference(tree2)


def get_wrf_distance(tree1, tree2):
    return tree1.robinson_foulds_distance(tree2)


def get_euc_distance(tree1, tree2):
    return tree1.euclidean_distance(tree2)


def convert_to_dendropy_trees(trees):
    taxa = dpy.TaxonSet()
    dpy_tree_list = [dpy.Tree.get_from_string(tree.newick, 'newick',
                     taxon_set=taxa) for tree in trees]
    return dpy_tree_list