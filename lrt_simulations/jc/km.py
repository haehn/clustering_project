#!/usr/bin/python
# -*- coding: utf-8 -*-
# Playing with affinity propagation from scikit-learn

import numpy as np
from sklearn.manifold import MDS
from sklearn.cluster import KMeans
import cPickle
import pylab as pl
from mpl_toolkits.mplot3d import Axes3D

# Load the sequence collection for the Rokas dataset
yg = cPickle.load(file('sim4_jc69_geodesic.pickle'))

# Geodesic distance matrix
distance_matrix = yg.get_distance_matrices().values()[0]

def double_centre(matrix):
    size = len(matrix)
    output = np.zeros( [size,size] )
    row_means = np.array([np.mean(row) for row in matrix])
    col_means = np.array([np.mean(col) for col in matrix.T])
    matrix_mean = np.mean(matrix)
    for i in range(size):
        for j in range(size):
            output[i,j] = matrix[i,j] - row_means[i] - col_means[j] + matrix_mean

    return output / -2.

# Embed the dm in 2 dimensions with multidimensional scaling
embed = MDS(n_components=100)
coords = embed.fit_transform(distance_matrix)
# coords = double_centre(distance_matrix)
# coords = distance_matrix
colors = {0:'r',1:'b',2:'g',3:'y'}

#KMeans
est = KMeans(n_clusters=4)
fig = pl.figure()
ax = fig.add_subplot(111, projection='3d')
est.fit(coords)
labels = est.labels_
print labels
for i in range(len(coords)):
  ax.scatter(coords[i:i+1, 0], coords[i:i+1, 1], coords[i:i+1, 2], c=colors[labels[i]])

ax.w_xaxis.set_ticklabels([])
ax.w_yaxis.set_ticklabels([])
ax.w_zaxis.set_ticklabels([])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

pl.show()