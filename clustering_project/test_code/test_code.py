#!/usr/bin/python
# -*- coding: utf-8 -*-

from sequence_collection import SequenceCollection
import numpy as np
import copy
import cPickle
np.set_printoptions(precision=2, linewidth=200)

"""
Protocol for running analysis:

1)  Create SequenceCollection object, e.g.:
    col = SequenceCollection( directory, datatype = 'protein' )
    This reads sequence alignments from an input directory
    (fasta format, phylip sequential, or phylip interleaved)
    Datatypes allowed are 'protein' and 'dna'

2)  Get phylogenetic trees for each gene with
    col.put_trees_parallel( program='treecollection', tmpdir = '/tmp' )
    Available programs - 'phyml', 'raxml', 'treecollection'

3)  Do hierarchical clustering on the trees, e.g.:
    col.put_partitions( metrics=['euc', 'geodesic', 'sym'], linkages=['single','complete'],
        nclasses=[4,5,6])
    available metrics are 'euc' for euclidean distance, 'geodesic', 'rf' for 
    weighted Robinson-Foulds (with branch lengths), and 'sym' for symmetric difference,
    i.e. standard, topology-only Robinson-Foulds

4)  Propagate the results through the data structure with:
    col.put_clusters()
    col.put_cluster_trees_parallel()
    
5)  Compare scores derived from clusters to random permutation of the original data
    either by making a copy of the SequenceCollection object, with clusters made up
    of the same number of genes with the same number of characters, or by randomising
    the alignments and performing hierarchical clustering on the randomised data

    if the former, do rand1 = col.make_randomised_copy
    if the latterm do rand2 = SequenceCollection(records=col.get_randomised_alignments(),
        datatype = 'protein')
"""


# indir = '/Users/kgori/git/kevin/yeast_data/MSA'
indir = '/Users/kgori/git/kevin/data/simulated_data/eight/MSA'

col = SequenceCollection(indir, datatype='protein')
ran = SequenceCollection(records=col.get_randomised_alignments(), datatype='protein')
col.put_trees_parallel()
ran.put_trees_parallel()
col.put_partitions(metrics=['euc','rf','sym'], linkages=['ward'], nclasses=[2,3,4,5,6,7,8,9,10])
ran.put_partitions(metrics=['euc','rf','sym'], linkages=['ward'], nclasses=[2,3,4,5,6,7,8,9,10])
col.put_clusters()
col.put_cluster_trees_parallel()
ran.put_clusters()
ran.put_cluster_trees_parallel()
rn2 = col.make_randomised_copy()

r1 = ran.get_clusters()
r2 = rn2.get_clusters()
cl = col.get_clusters()

for k in sorted(cl):
    print k
    print 'Clustering from true data:                     ', cl[k].score
    print 'Clustering from randomised data:               ', r1[k].score
    print 'Clustering from true data + randomising result:', r2[k].score
    print
  
cPickle.dump(col, file('col.pickle','w'))
cPickle.dump(ran, file('ran.pickle','w'))
cPickle.dump(rn2, file('rn2.pickle','w'))

"""
YEAST RESULT
('euc', 'ward', 2)
Clustering from true data:                      564.248
Clustering from randomised data:                503.7414
Clustering from true data + randomising result: 505.216

('euc', 'ward', 3)
Clustering from true data:                      559.109
Clustering from randomised data:                502.2794
Clustering from true data + randomising result: 505.216

('euc', 'ward', 4)
Clustering from true data:                      559.1094
Clustering from randomised data:                502.2794
Clustering from true data + randomising result: 504.764

('euc', 'ward', 5)
Clustering from true data:                      558.8582
Clustering from randomised data:                501.863
Clustering from true data + randomising result: 504.7643

('euc', 'ward', 6)
Clustering from true data:                      552.39
Clustering from randomised data:                497.3842
Clustering from true data + randomising result: 504.7646

('euc', 'ward', 7)
Clustering from true data:                      552.3895
Clustering from randomised data:                497.2594
Clustering from true data + randomising result: 504.1181

('euc', 'ward', 8)
Clustering from true data:                      551.3077
Clustering from randomised data:                479.8303
Clustering from true data + randomising result: 501.5788

('geodesic', 'ward', 2)
Clustering from true data:                      566.453
Clustering from randomised data:                505.511
Clustering from true data + randomising result: 507.962

('geodesic', 'ward', 3)
Clustering from true data:                      562.8588
Clustering from randomised data:                503.738
Clustering from true data + randomising result: 507.9623

('geodesic', 'ward', 4)
Clustering from true data:                      560.6918
Clustering from randomised data:                485.9653
Clustering from true data + randomising result: 507.1153

('geodesic', 'ward', 5)
Clustering from true data:                      560.6914
Clustering from randomised data:                485.9649
Clustering from true data + randomising result: 506.0329

('geodesic', 'ward', 6)
Clustering from true data:                      554.2232
Clustering from randomised data:                481.2518
Clustering from true data + randomising result: 506.0332

('geodesic', 'ward', 7)
Clustering from true data:                      552.7076
Clustering from randomised data:                479.7468
Clustering from true data + randomising result: 505.8623

('geodesic', 'ward', 8)
Clustering from true data:                      552.2001
Clustering from randomised data:                479.25047
Clustering from true data + randomising result: 502.6491

('sym', 'ward', 2)
Clustering from true data:                      528.312
Clustering from randomised data:                466.111
Clustering from true data + randomising result: 507.962

('sym', 'ward', 3)
Clustering from true data:                      494.925
Clustering from randomised data:                452.9
Clustering from true data + randomising result: 507.963

('sym', 'ward', 4)
Clustering from true data:                      476.5304
Clustering from randomised data:                445.2736
Clustering from true data + randomising result: 507.9634

('sym', 'ward', 5)
Clustering from true data:                      474.4247
Clustering from randomised data:                443.04789
Clustering from true data + randomising result: 506.1102

('sym', 'ward', 6)
Clustering from true data:                      471.6305
Clustering from randomised data:                443.04829
Clustering from true data + randomising result: 504.1927

('sym', 'ward', 7)
Clustering from true data:                      467.0354
Clustering from randomised data:                440.96919
Clustering from true data + randomising result: 504.1928

('sym', 'ward', 8)
Clustering from true data:                      465.6466
Clustering from randomised data:                440.36436
Clustering from true data + randomising result: 503.28341

SMALL
('euc', 'ward', 2)
Clustering from true data:                      6475.8809
Clustering from randomised data:                1861.68
Clustering from true data + randomising result: 1920.436

('euc', 'ward', 3)
Clustering from true data:                      2731.1542
Clustering from randomised data:                1842.685
Clustering from true data + randomising result: 1914.138

('euc', 'ward', 4)
Clustering from true data:                      170.0565
Clustering from randomised data:                1810.53
Clustering from true data + randomising result: 1910.846

('euc', 'ward', 5)
Clustering from true data:                      167.93243
Clustering from randomised data:                1809.101
Clustering from true data + randomising result: 1907.7012

('euc', 'ward', 6)
Clustering from true data:                      167.93233
Clustering from randomised data:                1800.9
Clustering from true data + randomising result: 1902.0062

('euc', 'ward', 7)
Clustering from true data:                      167.93233
Clustering from randomised data:                1800.9
Clustering from true data + randomising result: 1901.7676

('euc', 'ward', 8)
Clustering from true data:                      167.22473
Clustering from randomised data:                1797.0312
Clustering from true data + randomising result: 1898.6548

('geodesic', 'ward', 2)
Clustering from true data:                      6475.8809
Clustering from randomised data:                1882.06
Clustering from true data + randomising result: 1920.436

('geodesic', 'ward', 3)
Clustering from true data:                      3682.3542
Clustering from randomised data:                1842.685
Clustering from true data + randomising result: 1914.138

('geodesic', 'ward', 4)
Clustering from true data:                      170.0565
Clustering from randomised data:                1810.53
Clustering from true data + randomising result: 1910.846

('geodesic', 'ward', 5)
Clustering from true data:                      167.93243
Clustering from randomised data:                1802.329
Clustering from true data + randomising result: 1907.7012

('geodesic', 'ward', 6)
Clustering from true data:                      167.93233
Clustering from randomised data:                1798.4602
Clustering from true data + randomising result: 1902.0062

('geodesic', 'ward', 7)
Clustering from true data:                      167.93233
Clustering from randomised data:                1797.2582
Clustering from true data + randomising result: 1901.7676

('geodesic', 'ward', 8)
Clustering from true data:                      165.55163
Clustering from randomised data:                1787.1609
Clustering from true data + randomising result: 1905.9766

('sym', 'ward', 2)
Clustering from true data:                      6475.8909
Clustering from randomised data:                1835.746
Clustering from true data + randomising result: 1920.436

('sym', 'ward', 3)
Clustering from true data:                      2730.8042
Clustering from randomised data:                1814.139
Clustering from true data + randomising result: 1914.138

('sym', 'ward', 4)
Clustering from true data:                      170.0565
Clustering from randomised data:                1812.305
Clustering from true data + randomising result: 1910.846

('sym', 'ward', 5)
Clustering from true data:                      168.59519
Clustering from randomised data:                1808.871
Clustering from true data + randomising result: 1908.525

('sym', 'ward', 6)
Clustering from true data:                      166.21449
Clustering from randomised data:                1800.165
Clustering from true data + randomising result: 1912.732

('sym', 'ward', 7)
Clustering from true data:                      164.22559
Clustering from randomised data:                1791.6549
Clustering from true data + randomising result: 1886.527

('sym', 'ward', 8)
Clustering from true data:                      163.38334
Clustering from randomised data:                1790.0077
Clustering from true data + randomising result: 1885.2419
"""