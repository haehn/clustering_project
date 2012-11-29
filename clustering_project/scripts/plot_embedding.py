#!/usr/bin/env python

from sequence_collection import SequenceCollection
import os
import argparse
import sys
import glob
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import cPickle
import rpy2.robjects as rob
import uuid
print 'imports complete'


def fpath(s):
    """
    Helper function used when passing filepath arguments with argparse module.
    Trims all '/' characters from the end of the path string.
    """

    while s.endswith('/'):
        s = s[:-1]
    return s


def centre_of_mass(coords):
    return np.mean(coords, axis=0)


def get_coords(dm):
    dbc = sc.Clustering.get_double_centre(dm)
    (eigvals, eigvecs, cve) = sc.Clustering.get_eigen(dbc,
            standardize=True)
    (coords, varexp) = sc.Clustering.get_coords_by_dimension(eigvals,
            eigvecs, cve, 3, normalise=False)
    return coords


desc = \
    '''Plots an embedding of the trees in a Principal Coordinate space,
and saves as pdf. Also plots summary trees as pdfs, and a distance matrix'''
parser = argparse.ArgumentParser(prog=sys.argv[0], description=desc,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-i', '--input_directory',
                    help='Directory to work on', type=fpath)

args = vars(parser.parse_args())
indir = args['input_directory']
gtp_path = \
    '/nfs/research2/goldman/kevin/git/kevin/clustering_project/class_files'
helper = \
    '/nfs/research2/goldman/kevin/git/kevin/clustering_project/class_files/DV_wrapper.drw'
uid = uuid.uuid1()
while os.path.isdir('/nfs/nobackup/research/goldman/kg/{0}'.format(uid)):
    uid = uuid.uuid1()
tmpdir = '/nfs/nobackup/research/goldman/kg/{0}'.format(uid)
os.mkdir(tmpdir)

# gtp_path = os.environ['GTP_PATH']
# helper = os.environ['DARWINHELPER']
# if 'TEMPORARY_DIRECTORY' in os.environ:
#     tmpdir = os.environ['TEMPORARY_DIRECTORY']
# else:
#     tmpdir = '/tmp'
#     print 'Defaulting the tmpdir'

print indir

rob.r('library(ape)')
rob.r('library(phangorn)')
print 'r libraries loaded'

alignments_dir = '{0}/dna_alignments'.format(indir)

if os.path.isfile('{0}/tmpPickle.pkl'.format(indir)):
    sc = cPickle.load(file('{0}/tmpPickle.pkl'.format(indir)))
else:

    sc = SequenceCollection(
        alignments_dir,
        file_format='fasta',
        datatype='dna',
        gtp_path=gtp_path,
        helper=helper,
        tmpdir=tmpdir,
        )

    sc.put_trees(program='bionj')
    sc.put_partitions('geo', 'spectral', 4)
    sc.concatenate_records()
    sc.put_cluster_trees(program='bionj')
    cPickle.dump(sc, open('{0}/tmpPickle.pkl'.format(indir), 'w'))
print 'SC object available'

# Plot the heatmap of the distance matrix

dm = sc.get_distance_matrices()['geo']
p = sc.partitions[sc.clusters_to_partitions[('geo', 'spectral', 4)]]
fig1 = dm.plot_heatmap(sort_partition=p.get_memberships(flatten=True))
fig1.savefig('{0}/spectral_distance_matrix.pdf'.format(indir))
print 'distance matrix plot complete'

# Plot the embedding

coords = get_coords(dm.matrix)
min_Z = min([z for x,y,z in coords])
p = np.array(p.partition_vector)
p1 = np.where(p == 1)
p2 = np.where(p == 2)
p3 = np.where(p == 3)
p4 = np.where(p == 4)

cluster_trees = dict((t.name, t) for t in sc.get_cluster_trees())

colors = 'bgrcmyk'
coldict = {
    'b': 'blue',
    'g': 'green',
    'r': 'red',
    'c': 'cyan',
    'm': 'magenta',
    'y': 'yellow',
    'k': 'black',
    }
fig2d = plt.figure()
fig3d = plt.figure()
ax2d = fig2d.add_subplot(111)
ax3d = fig3d.add_subplot(111, projection='3d')

trees = {}

for (pos, partition) in enumerate((p1, p2, p3, p4)):
    hook = '-'.join([str(x) for x in partition[0]])
    tree = cluster_trees[hook].newick
    rob.r('t<-read.tree(text="' + tree + '")')
    rob.r('m <- midpoint(t)')
    fname = '{0}/cluster{1}_{2}.pdf'.format(indir, pos + 1,
            coldict[colors[pos]])
    rob.r('pdf("' + fname + '")')
    rob.r('plot(m, edge.width = 3, cex = 1.5, label.offset = 0.02)')
    rob.r('dev.off()')

    for i in partition[0]:
        ax2d.scatter(color=colors[pos], *(coords[i])[:2])
        ax3d.scatter(color=colors[pos], *coords[i])
        ax3d.plot(
            [coords[i][0], coords[i][0]],
            [coords[i][1], coords[i][1]],
            [min_Z, coords[i][2]],
            color=colors[pos],
            linewidth=0.2,
            alpha=0.7,
            )
    com = centre_of_mass(coords[partition])
    ax2d.scatter(color='k', marker='x', s=2, *com[:2])
    ax3d.scatter(color='k', marker='x', s=2, *com)

ax2d.set_xlabel('PCo1')
ax2d.set_ylabel('PCo2')
ax2d.set_title('Trees embedded in dimension-reduced space')
ax3d.set_xlabel('PCo1')
ax3d.set_ylabel('PCo2')
ax3d.set_zlabel('PCo3')
ax3d.set_title('Trees embedded in dimension-reduced space')
plt.show()
fig2d.savefig('{0}/embedding2d.pdf'.format(indir))
fig3d.savefig('{0}/embedding3d.pdf'.format(indir))

# tests

try:
    for f in ['{0}/embedding2d.pdf'.format(indir),
              '{0}/embedding3d.pdf'.format(indir),
              '{0}/spectral_distance_matrix.pdf'.format(indir)]:
        assert os.path.isfile(f)
except AssertionError:
    print 'Not all pdf files were created'
    sys.exit(1)
# os.system('rm {0}/tmpPickle.pkl'.format(indir))
os.system('rm -r {0}'.format(tmpdir))
