#!/usr/bin/env python
from ete2 import Tree, TreeFace, TreeStyle, NodeStyle, faces, AttrFace
import glob, sys
from handleArgs import handleArgs

def layout(node):
    if node.is_leaf():
        N = AttrFace("name", fsize=14, fgcolor="black")
        #faces.add_face_to_node(N, node, 0)
        small_tree = small_d[node.name]
        T = TreeFace(small_tree,small_ts)
        T.opacity = 0.8
        # And place as a float face over the tree
        faces.add_face_to_node(T, node, 1, position="aligned")

args = handleArgs(sys.argv,help='''
Provide the input directory and optional rooting for Trees
-dir  =  Input dir
-root =  root
''')

if "-dir" not in args:
    INPUT_DIR = '.'
elif not args['-dir']:
    print "No input directory specified, trying '.' ..."
    INPUT_DIR = ','
else: INPUT_DIR = args['-dir']

if "-root" not in args:
    ROOT = False
elif not args["-root"]:
    ROOT = False
else: ROOT = args["-root"]

tree_files = glob.glob( "{0}/trees/besttrees/*.nwk".format(INPUT_DIR))
small = [Tree( open(x).read() ) for x in tree_files]
for x in small:
    x.set_outgroup(ROOT)
    x.ladderize()
names = [x[1+x.rindex("/"):x.rindex(".")] for x in tree_files]
big = Tree( open("{0}/clusters/cluster_dendrogram.nwk".format(INPUT_DIR)).read() )
small_d = dict(zip(names,small))
small_ts = TreeStyle()
small_ts.tree_width = 200
small_ts.show_leaf_name = True
# small_ts.rotation = -90
big_ts = TreeStyle()
big_ts.show_leaf_name = True#False
big_ts.layout_fn = layout
# big_ts.rotation = 90
big_ts.mode = 'r'


big.show(tree_style=big_ts)
