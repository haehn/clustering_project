#!/usr/bin/python
# -*- coding: utf-8 -*-
# Playing with affinity propagation from scikit-learn

import numpy as np
from sklearn.manifold import MDS
from sklearn.cluster import AffinityPropagation
import cPickle
import pylab as pl
from itertools import cycle
from mpl_toolkits.mplot3d import Axes3D

# Load the sequence collection for the Rokas dataset
yg = cPickle.load(file('yeast_jc69_geodesic.pickle'))

# Geodesic distance matrix
distance_matrix = yg.get_distance_matrices().values()[0]

# Embed the dm in 2 dimensions with multidimensional scaling
embed = MDS(n_components=3)
coords = embed.fit_transform(distance_matrix)

# Recalculate dm from embedded coords
coords_norms = np.sum(coords ** 2, axis=1)
S = - coords_norms[:, np.newaxis] - coords_norms[np.newaxis, :] + 2 * np.dot(coords, coords.T)

# Affinity propagation
p = 5 * np.median(S)
af = AffinityPropagation().fit(S)
cluster_centers_indices = af.cluster_centers_indices_
labels = af.labels_
print 'Labels = %s' % labels

n_clusters_ = len(cluster_centers_indices)

print 'Estimated number of clusters: %d' % n_clusters_

# Plot

pl.close('all')
pl.figure(1)
pl.clf()

colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
for k, col in zip(range(n_clusters_), colors):
    class_members = labels == k
    cluster_center = coords[cluster_centers_indices[k]]
    pl.plot(coords[class_members, 0], coords[class_members, 1], col + '.')
    pl.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
            markeredgecolor='k', markersize=14)
    for x in coords[class_members]:
        pl.plot([cluster_center[0], x[0]], [cluster_center[1], x[1]], col)

pl.title('Estimated number of clusters: %d' % n_clusters_)
pl.show()

