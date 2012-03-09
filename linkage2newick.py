#!/usr/bin/env python

import numpy as np
l = np.array([[  6.,     9.,     0.,     2.  ],
 [  5.,     7.,     0.,     2.  ],
 [  1.,     4.,     0.,     2.  ],
 [ 10.,    11.,     0.17,   4.  ],
 [  3.,    12.,     0.17,   3.  ],
 [  8.,    13.,     0.33,   5.  ],
 [  2.,    14.,     0.33,   4.  ],
 [  0.,    16.,     0.5,    5.  ],
 [ 15.,    17.,     1.,    10.  ]])

l2 = np.array([[  1.,     4.,     0.11,   2.  ],
 [  7.,     8.,     0.14,   2.  ],
 [  6.,     9.,     0.14,   2.  ],
 [  0.,    10.,     0.17,   3.  ],
 [  5.,    12.,     0.18,   3.  ],
 [  2.,     3.,     0.18,   2.  ],
 [ 13.,    15.,     0.19,   5.  ],
 [ 11.,    14.,     0.3,    5.  ],
 [ 16.,    17.,     1.,    10.  ]])

l3 = np.array([[ 12.,          19.,           0.,           1.        ],
[ 15.,          21.,           0.,           3.        ],
[ 18.,          22.,           0.,           4.        ],
[  3.,          16.,           0.,           7.        ],
[  8.,          23.,           0.,           6.        ],
[  5.,          27.,           0.,           6.        ],
[  1.,          28.,           0.,           7.        ],
[  0.,          21.,           0.,           2.        ],
[  5.,          29.,           0.18350472,   2.        ],
[  2.,          10.,           0.18350472,   3.        ],
[ 47.,          30.,           0.29289577,   9.        ],
[ 13.,          28.,           0.29289577,  13.        ],
[ 73.,          32.,           0.29289577,  18.        ],
[ 26.,          12.,           0.42264521,   5.        ],
[  5.,          33.,           0.42264521,  12.        ],
[ 14.,          35.,           0.42264521,  12.        ],
[ 19.,          35.,           0.42264521,  18.        ],
[  4.,          20.,           0.31174826,   3.        ],
[ 34.,          21.,           0.5,         19.        ],
[ 38.,          29.,           0.31174826,  21.        ]])

g3 = ["l{0}".format(i+1) for i in range(len(l3)+1)]

g = ["gene001","gene002","gene003","gene004","gene005","gene006","gene007","gene008","gene009","gene010"]

def linkage2newick(gene_list,linkage):
    import re
    d = {}
    num_genes = len(gene_list)
    for i in range(len(linkage)+len(gene_list)):
        if i < num_genes:
            d[i]=[gene_list[i],0]
        else:
            d[i] = [ (int(linkage[i-num_genes][0]),int(linkage[i-num_genes][1])),linkage[i-num_genes][2] ]

    s = str(i)+":0.0;"

    def quick_replace(s,d,i):
        # Replace node i with the subtree it represents
        subtree = d[i][0]
        if type(subtree)==type( (1,) ): #check tuple-ness before unpacking
            left, right = subtree
            left_brlen = d[i][1] - d[left][1] #calculate branch lengths for the unpacked subtrees from the tuple
            right_brlen = d[i][1] - d[right][1]
            out = re.sub( "(?<![\w\d\:.]){0}(?!\w\d)".format(i), "({0}:{1}, {2}:{3})".format(left,left_brlen,right,right_brlen), s ) #replace node with subtree
            return out
        else: # if tuple-ness fails then this is a leaf node, so replace with leaf name from genelist
            return re.sub( "(?<![\w\d\:.]){0}(?!\w\d)".format(i), "{0}".format(d[i][0]), s )

    # i is set to the number of the final node in the tree - this represents the whole tree contracted to a single node
    # iteratively replacing the node with the next-most contracted subtree will result in a fully expanded tree
    while i >= 0: # take advantage of i's current value at len(linkage)+len(gene_list), i.e. max
        s = quick_replace(s,d,i)
        i -= 1 # repeatedly call quick_replace function to rename node

    return s


linkage2newick(g,l)
linkage2newick(g3,l3)
