#!/usr/bin/env python

import dendropy as dpy


# DENDROPY UTILS

def convert_to_dendropy_tree(tree, taxon_set=None):
    if taxon_set is None:
        return dpy.Tree.get_from_string(tree.newick, 'newick')
    return dpy.Tree.get_from_string(tree.newick, 'newick',
                                    taxon_set=taxon_set)


def convert_to_dendropy_trees(trees):
    taxa = dpy.TaxonSet()
    return [convert_to_dendropy_tree(tree, taxa) for tree in trees]


def check_rooted(newick):
    t = dpy.Tree.get_from_string(newick, 'newick')
    root_degree = len(t.seed_node.child_nodes())
    return root_degree == 2


def deroot(newick):
    t = dpy.Tree.get_from_string(newick, 'newick')
    t.deroot()
    return t.as_newick_string() + ';\n'


def get_rf_distance(tree1, tree2):
    return tree1.symmetric_difference(tree2)


def get_wrf_distance(tree1, tree2):
    return tree1.robinson_foulds_distance(tree2)


def get_euc_distance(tree1, tree2):
    return tree1.euclidean_distance(tree2)
