import sys
import glob
import re
import os
import scipy.stats as ss
from pylab import *

sort_key = lambda item: tuple((int(num) if num else alpha) for (num,
                              alpha) in re.findall(r'(\d+)|(\D+)',
                              item))

basedir = sys.argv[1]
while basedir.endswith('/'):
    basedir = basedir[:-1]

simdirs = [g for g in glob.glob('{0}/*'.format(basedir))
           if os.path.isdir(g)]
simdirs.sort(key=sort_key)


def parse_results(d):
    symbase = None
    syminter = None
    symintra = None
    symover = None
    eucbase = None
    eucinter = None
    eucintra = None
    eucover = None
    symbase = None
    syminter = None
    symintra = None
    symover = None
    truescore = None
    symrecord = False
    eucrecord = False
    symrecord = False

    getval = lambda x: x.rstrip().split()[-1]
    with open('{0}/treedistances.txt'.format(d)) as file:
        for line in file:
            if line.startswith('geodesic'):
                geobase = float(getval(line))
            elif line.startswith('euclidean'):
                eucbase = float(getval(line))
            elif line.startswith('RF') \
                and not line.startswith('RF class'):
                symbase = float(getval(line))
            elif line.startswith('Geodesic class stats'):
                georecord = True
                eucrecord = False
                symrecord = False
            elif line.startswith('Euc class stats'):
                georecord = False
                eucrecord = True
                symrecord = False
            elif line.startswith('RF class stats'):
                georecord = False
                eucrecord = False
                symrecord = True
            elif line.startswith('wRF'):
                georecord = False                
                eucrecord = False
                symrecord = False
            elif line.startswith('inter_mean'):
                if georecord:
                    geointer = float(getval(line))
                elif eucrecord:
                    eucinter = float(getval(line))
                elif symrecord:
                    syminter = float(getval(line))
            elif line.startswith('intra_mean'):
                if georecord:
                    geointra = float(getval(line))
                elif eucrecord:
                    eucintra = float(getval(line))
                elif symrecord:
                    symintra = float(getval(line))
            elif line.startswith('overall_mean'):
                if georecord:
                    geoover = float(getval(line))
                elif eucrecord:
                    eucover = float(getval(line))
                elif symrecord:
                    symover = float(getval(line))
            elif line.startswith('True score'):
                truescore = float(getval(line))

    retdic = {
        'geo': {
            'inter': geointer,
            'intra': geointra,
            'overall': geoover,
            'base': geobase,
            },
        'euc': {
            'inter': eucinter,
            'intra': eucintra,
            'overall': eucover,
            'base': eucbase,
            },
        'sym': {
            'inter': syminter,
            'intra': symintra,
            'overall': symover,
            'base': symbase,
            },
        'truescore': truescore,
        }

    files = glob.glob('{0}/euc*'.format(d)) \
        + glob.glob('{0}/geo*'.format(d)) \
        + glob.glob('{0}/sym*'.format(d))
    for f in files:
        with open(f) as file:
            for line in file:
                if line.startswith('Distance'):
                    k = getval(line)
                if line.startswith('Method'):
                    method = getval(line)
                if line.startswith('Score'):
                    score = float(getval(line))
                if line.startswith('Varinf'):
                    vi = float(getval(line))
        retdic[k][method] = {'score': score, 'vi': vi}

    return retdic


collect = []

for directory in simdirs:
    print 'parsing {0}'.format(directory)
    collect.append(parse_results(directory))

geocomplete=[]
geoaverage=[]
geoward=[]
geosingle=[]
geokmedoids=[]
geospectral=[]
geoMDS=[]

euccomplete=[]
eucaverage=[]
eucward=[]
eucsingle=[]
euckmedoids=[]
eucspectral=[]
eucMDS=[]

symcomplete=[]
symaverage=[]
symward=[]
symsingle=[]
symkmedoids=[]
symspectral=[]
symMDS=[]

for record in collect:
    geocomplete.append((record['geo']['inter']-record['geo']['intra'],(record['geo']['complete']['score']-record['truescore'])/record['truescore']))
    geoaverage.append((record['geo']['inter']-record['geo']['intra'],(record['geo']['average']['score']-record['truescore'])/record['truescore']))
    geoward.append((record['geo']['inter']-record['geo']['intra'],(record['geo']['ward']['score']-record['truescore'])/record['truescore']))
    geosingle.append((record['geo']['inter']-record['geo']['intra'],(record['geo']['single']['score']-record['truescore'])/record['truescore']))
    geoMDS.append((record['geo']['inter']-record['geo']['intra'],(record['geo']['MDS']['score']-record['truescore'])/record['truescore']))
    geokmedoids.append((record['geo']['inter']-record['geo']['intra'],(record['geo']['kmedoids']['score']-record['truescore'])/record['truescore']))
    geospectral.append((record['geo']['inter']-record['geo']['intra'],(record['geo']['spectral']['score']-record['truescore'])/record['truescore']))

    euccomplete.append((record['euc']['inter']-record['euc']['intra'],(record['euc']['complete']['score']-record['truescore'])/record['truescore']))
    eucaverage.append((record['euc']['inter']-record['euc']['intra'],(record['euc']['average']['score']-record['truescore'])/record['truescore']))
    eucward.append((record['euc']['inter']-record['euc']['intra'],(record['euc']['ward']['score']-record['truescore'])/record['truescore']))
    eucsingle.append((record['euc']['inter']-record['euc']['intra'],(record['euc']['single']['score']-record['truescore'])/record['truescore']))
    eucMDS.append((record['euc']['inter']-record['euc']['intra'],(record['euc']['MDS']['score']-record['truescore'])/record['truescore']))
    euckmedoids.append((record['euc']['inter']-record['euc']['intra'],(record['euc']['kmedoids']['score']-record['truescore'])/record['truescore']))
    eucspectral.append((record['euc']['inter']-record['euc']['intra'],(record['euc']['spectral']['score']-record['truescore'])/record['truescore']))

    symcomplete.append((record['sym']['inter']-record['sym']['intra'],(record['sym']['complete']['score']-record['truescore'])/record['truescore']))
    symaverage.append((record['sym']['inter']-record['sym']['intra'],(record['sym']['average']['score']-record['truescore'])/record['truescore']))
    symward.append((record['sym']['inter']-record['sym']['intra'],(record['sym']['ward']['score']-record['truescore'])/record['truescore']))
    symsingle.append((record['sym']['inter']-record['sym']['intra'],(record['sym']['single']['score']-record['truescore'])/record['truescore']))
    symMDS.append((record['sym']['inter']-record['sym']['intra'],(record['sym']['MDS']['score']-record['truescore'])/record['truescore']))
    symkmedoids.append((record['sym']['inter']-record['sym']['intra'],(record['sym']['kmedoids']['score']-record['truescore'])/record['truescore']))
    symspectral.append((record['sym']['inter']-record['sym']['intra'],(record['sym']['spectral']['score']-record['truescore'])/record['truescore']))




