#!/usr/bin/env python

import sys, glob, os, re
sort_key = lambda item: tuple((int(num) if num else alpha) for (num,
                              alpha) in re.findall(r'(\d+)|(\D+)',
                              item))

basedir = os.path.abspath('/homes/kgori/storage/regime2')

#nnidirs = glob.glob('{0}/nni*'.format(basedir))
nnidirs = ['{0}/nni{1}'.format(basedir,i) for i in range(2,6)]
simdirs = []
for nnidir in nnidirs:
    simdirs += glob.glob('{0}/sim*'.format(nnidir))

simdirs.sort(key=sort_key)
#print basedir
#print nnidirs
#print simdirs
warnings = []
for d in simdirs:
    print 'Checking {0}'.format(d)
    if not os.path.isfile('{0}/phyml_clustering/true_cluster1.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/true_cluster2.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/true_cluster3.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/true_cluster4.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/eucMDS4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/eucMDS4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/eucMDS4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/eucMDS4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/eucave4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/eucave4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/eucave4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/eucave4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/euccom4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/euccom4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/euccom4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/euccom4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/euckme4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/euckme4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/euckme4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/euckme4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/eucsin4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/eucsin4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/eucsin4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/eucsin4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/eucspe4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/eucspe4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/eucspe4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/eucspe4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/eucwar4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/eucwar4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/eucwar4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/eucwar4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/geoMDS4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/geoMDS4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/geoMDS4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/geoMDS4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/geoave4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/geoave4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/geoave4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/geoave4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/geocom4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/geocom4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/geocom4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/geocom4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/geokme4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/geokme4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/geokme4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/geokme4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/geosin4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/geosin4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/geosin4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/geosin4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/geospe4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/geospe4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/geospe4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/geospe4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/geowar4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/geowar4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/geowar4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/geowar4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/symMDS4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/symMDS4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/symMDS4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/symMDS4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/symave4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/symave4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/symave4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/symave4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/symcom4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/symcom4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/symcom4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/symcom4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/symkme4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/symkme4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/symkme4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/symkme4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/symsin4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/symsin4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/symsin4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/symsin4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/symspe4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/symspe4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/symspe4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/symspe4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/symwar4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/symwar4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/symwar4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/phyml_clustering/symwar4_4.phy'.format(d)):

        warnings.append( (d, 'phyml_clustering') )
        #print d, 'phyml_clustering'
    elif not os.path.isfile('{0}/bionj_clustering/true_cluster1.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/true_cluster2.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/true_cluster3.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/true_cluster4.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/eucMDS4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/eucMDS4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/eucMDS4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/eucMDS4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/eucave4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/eucave4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/eucave4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/eucave4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/euccom4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/euccom4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/euccom4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/euccom4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/euckme4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/euckme4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/euckme4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/euckme4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/eucsin4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/eucsin4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/eucsin4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/eucsin4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/eucspe4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/eucspe4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/eucspe4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/eucspe4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/eucwar4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/eucwar4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/eucwar4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/eucwar4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/geoMDS4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/geoMDS4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/geoMDS4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/geoMDS4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/geoave4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/geoave4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/geoave4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/geoave4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/geocom4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/geocom4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/geocom4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/geocom4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/geokme4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/geokme4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/geokme4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/geokme4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/geosin4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/geosin4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/geosin4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/geosin4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/geospe4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/geospe4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/geospe4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/geospe4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/geowar4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/geowar4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/geowar4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/geowar4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/symMDS4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/symMDS4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/symMDS4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/symMDS4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/symave4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/symave4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/symave4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/symave4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/symcom4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/symcom4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/symcom4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/symcom4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/symkme4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/symkme4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/symkme4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/symkme4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/symsin4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/symsin4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/symsin4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/symsin4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/symspe4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/symspe4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/symspe4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/symspe4_4.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/symwar4_1.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/symwar4_2.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/symwar4_3.phy'.format(d)) \
        or not os.path.isfile('{0}/bionj_clustering/symwar4_4.phy'.format(d)):
        warnings.append( (d, 'bionj_clustering') ) 
        #print d, 'bionj_clustering'

for x in warnings:
    #pass 
    os.system('bsub -o /dev/null bash /homes/kgori/research/clustering_project/tempdir_wrapper.sh python doclustering.py -d {0}'.format(x[0]))
for x in warnings: print x
print len(warnings)
