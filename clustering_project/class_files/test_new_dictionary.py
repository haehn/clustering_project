#!/usr/bin/env python

print 'importing...'
from sequence_collection import SequenceCollection
from distance_matrix import DistanceMatrix
from partition import Partition
from clustering import Clustering
import sys
import numpy as np
from pylab import *
print 'done.'


def print_dict(d):
    for k in sorted(d):
        print d


np.set_printoptions(linewidth=200, precision=3)
sc = SequenceCollection(
    '/Users/kgori/scratch/chk/aa_alignments/',
    get_distances=False,
    file_format='phylip',
    helper='/Users/kgori/git/kevin/clustering_project/class_files/DV_wrapper.drw'
        ,
    parallel_load=True,
    gtp_path='/Users/kgori/git/kevin/clustering_project/class_files/',
    tmpdir='/tmp',
    datatype='protein',
    )

sc_yeast = SequenceCollection(
    '/Users/kgori/scratch/yeast_MSA',
    get_distances=False,
    file_format='phylip',
    helper='/Users/kgori/git/kevin/clustering_project/class_files/DV_wrapper.drw'
        ,
    parallel_load=True,
    gtp_path='/Users/kgori/git/kevin/clustering_project/class_files/',
    tmpdir='/tmp/yeast',
    datatype='dna',
    )

sc.put_trees(program='bionj', overwrite=False)
sc.put_partitions(['geo', 'rf'], ['spectral','NJW','ShiMalik'], range(1, 41))
# sc.put_partitions(['geo', 'rf'], ['spectral', 'MDS', 'kmedoids','single','complete','average','ward','NJW'], range(1, 41))
# sc.put_partitions(['geo', 'rf'], ['spectral-prune'], range(1, 41), recalculate=True)

sc_yeast.put_trees(program='bionj', overwrite=False)
sc_yeast.put_partitions(['geo', 'rf'], ['spectral', 'MDS'], range(1, 50))

# yr = sc_yeast.make_randomised_copy(get_distances=False)
# yr.put_trees(program='bionj', overwrite=False)
# yr.put_partitions(['geo', 'rf'], ['spectral', 'MDS'], range(1, 50))


sc.concatenate_records()
sc.put_cluster_trees_parallel(program='bionj', overwrite=False)

sc_yeast.concatenate_records()
sc_yeast.put_cluster_trees(program='bionj', overwrite=False)

# yr.concatenate_records()
# yr.put_cluster_trees(program='bionj', overwrite=False)
def extract(sc, method, distance):
    BICdic = {}
    likdic = {}
    print '''


    ::: :::


    '''
    for k in sorted(sc.clusters_to_partitions, key=lambda x: (x[0], x[2],
                    x[1])):
        if method in k and distance in k:
            score = sc.partitions[sc.clusters_to_partitions[k]].score
            nspec = sc.records[0].length
            sample_size = sc.length
            free_params = k[2] * 3 * (nspec - 2)
            BIC = -2 * score + free_params * np.log(sample_size)
            print k, '\t', score, '\t', BIC
            BICdic[k[2]] = BIC
            likdic[k[2]] = score

    return {'BIC':BICdic, 'lnL':likdic}

def plottables(d):
    z1 = d.items()
    x1 = [a for (a, b) in z1]
    y1 = [b for (a, b) in z1]
    y1_alt = [y1[i]-y1[i-1] for i in range(1,len(y1))]

    return x1, y1, y1_alt


# dMDS=extract(sc, 'MDS', 'geo')
# dMDSrf=extract(sc, 'MDS', 'rf')
dspectral=extract(sc, 'spectral', 'geo')
dspectralrf=extract(sc, 'spectral', 'rf')
# dkmedoids=extract(sc, 'kmedoids', 'geo')
# dkmedoidsrf=extract(sc, 'kmedoids', 'rf')
# dward=extract(sc, 'ward', 'geo')
# dwardrf=extract(sc, 'ward', 'rf')
# dsingle=extract(sc, 'single', 'geo')
# dsinglerf=extract(sc, 'single', 'rf')
# dcomplete=extract(sc, 'complete', 'geo')
# dcompleterf=extract(sc, 'complete', 'rf')
# daverage=extract(sc, 'average', 'geo')
# daveragerf=extract(sc, 'average', 'rf')
dNJW=extract(sc, 'NJW', 'geo')
dNJWrf=extract(sc, 'NJW', 'rf')
dShiMalik=extract(sc, 'ShiMalik', 'geo')
dShiMalikrf=extract(sc, 'ShiMalik', 'rf')

(xspectral,yspectral, y2spectral ) = plottables(dspectral['lnL'])
(xspectralrf,yspectralrf, y2spectralrf ) = plottables(dspectral['lnL'])
# (xMDS,yMDS, y2MDS ) = plottables(dMDS['lnL'])
# (xMDSrf,yMDSrf, y2MDSrf ) = plottables(dMDSrf['lnL'])
# (xkmedoids,ykmedoids, y2kmedoids ) = plottables(dkmedoids['lnL'])
# (xkmedoidsrf,ykmedoidsrf, y2kmedoidsrf ) = plottables(dkmedoidsrf['lnL'])
# (xward,yward, y2ward ) = plottables(dward['lnL'])
# (xwardrf,ywardrf, y2wardrf ) = plottables(dwardrf['lnL'])
# (xsingle,ysingle, y2single ) = plottables(dsingle['lnL'])
# (xsinglerf,ysinglerf, y2singlerf ) = plottables(dsinglerf['lnL'])
# (xcomplete,ycomplete, y2complete ) = plottables(dcomplete['lnL'])
# (xcompleterf,ycompleterf, y2completerf ) = plottables(dcompleterf['lnL'])
# (xaverage,yaverage, y2average ) = plottables(daverage['lnL'])
# (xaveragerf,yaveragerf, y2averagerf ) = plottables(daveragerf['lnL'])
(xNJW,yNJW, y2NJW ) = plottables(dNJW['lnL'])
(xNJWrf,yNJWrf, y2NJWrf ) = plottables(dNJWrf['lnL'])
(xShiMalik,yShiMalik, y2ShiMalik ) = plottables(dShiMalik['lnL'])
(xShiMalikrf,yShiMalikrf, y2ShiMalikrf ) = plottables(dShiMalikrf['lnL'])




# d2 = {}
# print '''


# :::YEAST:::


# '''
# for k in sorted(sc_yeast.clusters_to_partitions, key=lambda x: (x[0],
#                 x[2], x[1])):
#     if 'MDS' in k and 'geo' in k:
#         score = \
#             sc_yeast.partitions[sc_yeast.clusters_to_partitions[k]].score
#         free_params = k[2] * 18
#         BIC = -2 * score + free_params * np.log(106)
#         print k, score, BIC
#         d2[k[2]] = BIC
# print '''


# :::RANDOMISED YEAST:::


# '''
# d3 = {}
# for k in sorted(yr.clusters_to_partitions, key=lambda x: (x[0],
#                 x[2], x[1])):
#     if 'MDS' in k and 'geo' in k:
#         score = \
#             yr.partitions[yr.clusters_to_partitions[k]].score
#         free_params = k[2] * 18
#         BIC = -2 * score + free_params * np.log(106)
#         print k, score, BIC
#         d3[k[2]] = BIC

# z1 = d.items()
# x1 = [a for (a, b) in z1]
# y1 = [b for (a, b) in z1]
# y1_alt = [y1[i]-y1[i-1] for i in range(1,len(y1))]

# z2 = d2.items()
# x2 = [a for (a, b) in z2]
# y2 = [b for (a, b) in z2]
# y2_alt = [y2[i]-y2[i-1] for i in range(1,len(y2))]

# z3 = d3.items()
# x3 = [a for (a, b) in z3]
# y3 = [b for (a, b) in z3]
# y3_alt = [y3[i]-y3[i-1] for i in range(1,len(y3))]

# plot(x1, y1, 'b-')
# show()

# plot(x2, y2, 'b-')
# show()

# plot(x3, y3, 'r-')
# show()

# plot(x2[1:], y2_alt, 'b-')
# plot(x2[1:], y3_alt, 'r-')
# axhline(0)
# show()
# # print sc.get_trees()
# # print trees
# # print trees == sc.get_trees()

# # sc.put_distance_matrices(['rf','wrf','geo','euc'])
# # print sc.distance_matrices

# # sc.put_partition('rf','MDS',4)
# # print sc.clusters_to_partitions
# # print sc.partitions
# # sys.exit()

# # C = Clustering()
# # dm = sc.distance_matrices['rf']
# # kmed = C.run_clustering(dm,4,'kmedoids')
# # spec = C.run_clustering(dm,4,'spectral')
# # mdsc = C.run_clustering(dm,4,'MDS')
# # sing = C.run_clustering(dm,4,'single')
# # comp = C.run_clustering(dm,4,'complete')
# # aver = C.run_clustering(dm,4,'average')
# # ward = C.run_clustering(dm,4,'ward')
# # dbsc = C.run_clustering(dm,4,'dbscan')

# # for x in [kmed,spec,mdsc,sing,comp,aver,ward,dbsc]:
# #     print x

# # print sys.getsizeof(sc)

# sc.clusters_to_partitions['true'] = (
#     1,
#     1,
#     1,
#     1,
#     1,
#     2,
#     2,
#     2,
#     2,
#     2,
#     2,
#     2,
#     2,
#     2,
#     3,
#     3,
#     3,
#     3,
#     3,
#     3,
#     3,
#     3,
#     3,
#     3,
#     3,
#     3,
#     4,
#     4,
#     4,
#     4,
#     4,
#     4,
#     4,
#     4,
#     4,
#     4,
#     4,
#     4,
#     4,
#     4,
#     )

# sc.partitions[(
#     1,
#     1,
#     1,
#     1,
#     1,
#     2,
#     2,
#     2,
#     2,
#     2,
#     2,
#     2,
#     2,
#     2,
#     3,
#     3,
#     3,
#     3,
#     3,
#     3,
#     3,
#     3,
#     3,
#     3,
#     3,
#     3,
#     4,
#     4,
#     4,
#     4,
#     4,
#     4,
#     4,
#     4,
#     4,
#     4,
#     4,
#     4,
#     4,
#     4,
#     )] = Partition((
#     1,
#     1,
#     1,
#     1,
#     1,
#     2,
#     2,
#     2,
#     2,
#     2,
#     2,
#     2,
#     2,
#     2,
#     3,
#     3,
#     3,
#     3,
#     3,
#     3,
#     3,
#     3,
#     3,
#     3,
#     3,
#     3,
#     4,
#     4,
#     4,
#     4,
#     4,
#     4,
#     4,
#     4,
#     4,
#     4,
#     4,
#     4,
#     4,
#     4,
#     ))

