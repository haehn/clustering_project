from collections import defaultdict
import cPickle
import glob
import re
import sys
from sequence_collection import SequenceCollection


def load_records(directory, extension='*.pickle'):
    return [cPickle.load(open(f)) for f in
            glob.iglob('{0}/{1}'.format(directory, extension))]

def make_target_dict(target_list):
    return dict( (rec.name, rec) for rec in target_list)

def order(l):
    """
    The clustering returned by the hcluster module gives 
    group membership without regard for numerical order 
    This function preserves the group membership, but sorts 
    the labelling into numerical order
    
    *Faster version without recursion*
    """

    list_length = len(l)

    d = defaultdict(list)
    for (i, element) in enumerate(l):
        d[element].append(i)

    l2 = [None]*list_length

    for (name, index_list) in enumerate(sorted(d.values(), key=min), start=1):
        for index in index_list:            
            l2[index] = name

    return l2

def rebuild_partitions(query_list, target_dict, metric, method):
    """
    The query list is the individual alignments record list,
    The target list is the clustered alignments record list
    """

    partition_list = [0] * len(query_list)
    partition_list_check_copy = [0] * len(query_list)
    # print 'partition_list=', partition_list
    # print 'partition_list_check_copy', partition_list_check_copy
    header = query_list[0].headers[0]
    header_check_copy = query_list[0].headers[-1]
    # print 'header=',header
    # print 'header_check_copy=', header_check_copy
    key = metric[:3] + method[:3] + '4_'
    # print 'key=',key

    for i in range(1,5):
        # print 'i=',i
        for j, rec in enumerate(query_list):
            # print 'j=',j
            query_seq = rec.mapping[header]
            query_seq_check_copy = rec.mapping[header_check_copy]
            if query_seq in target_dict[key+str(i)].mapping[header]:
                partition_list[j] = i
            if query_seq_check_copy in target_dict[key+str(i)].mapping[header_check_copy]:
                partition_list_check_copy[j] = i

    try:
        assert partition_list == partition_list_check_copy
    except AssertionError:
        print 'Can\'t get partitions to match'
        return -1
    return partition_list

sort_key = lambda item: tuple((int(num) if num else alpha) for (num,alpha) in re.findall(r'(\d+)|(\D+)', item.name))

print 'loading recs...'
recs = sorted(load_records(sys.argv[1], sys.argv[2]), key=sort_key)
print 'loading crecs...'
crecs = load_records(sys.argv[3])
d = make_target_dict(crecs)
print 'done.'

t = [1]*15+[2]*15+[3]*15+[4]*15
c = SequenceCollection()
c.records = recs

for metric in ['sym', 'euc', 'geo']:
    c.put_distance_matrices(metric, gtp_path = '/homes/kgori/research/clustering_project/class_files')
    for method in ['single', 'complete', 'average', 'ward', 'MDS', 'spectral', 'kmedoids']:
        p = order(rebuild_partitions(recs, d, metric, method))
        c.clustering.partitions[(metric, method, 4)] = p
c.clustering.partitions['true'] = t
c.put_clusters()
print '(Done).'

for rec in c.get_cluster_records():
    try:
        rec.tree = d[rec.name].tree
    except KeyError:
        new_key = 'true_cluster{0}'.format(rec.name[-1])
        rec.tree = d[new_key].tree
c.get_clusters()['true'].update_score()
for metric in ['sym', 'euc', 'geo']:       
    for method in ['single', 'complete', 'average', 'ward', 'MDS', 'spectral', 'kmedoids']:
        p = c.clustering.partitions[(metric, method, 4)]
        res = c.get_clusters()[(metric, method, 4)]
        res.update_score()
        print '{0} {1: <4}{2} {3:.2f} {4:.0f}'.format(metric, method[:3], p, c.clustering.variation_of_information(p,t), res.score-c.get_clusters()['true'].score)

