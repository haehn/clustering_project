import glob
from pyplot import *
import numpy as np

fi = glob.glob('res4/*')
p = lambda x: x[x.rindex('/')+1:x.rindex('.')].split('_')[1:]
res4 = {}


for f in fi:
    metric, linkage = p(f)
    if not metric in res4:
        res4[metric] = {}
    if not linkage in res4:
        res4[linkage] = {}
    if not linkage in res4[metric]:
        res4[metric][linkage] = {'VI':[],'SoS':[]}
    if not metric in res4[linkage]:
        res4[linkage][metric] = {'VI':[],'SoS':[]}
    open_file = open(f)
    for line in open_file:
        VI, SoS = line.rstrip().split()
        VI = float(VI)
        SoS = float(SoS)
        res4[linkage][metric]['VI'].append(VI)
        res4[linkage][metric]['SoS'].append(SoS)
        res4[metric][linkage]['VI'].append(VI)
        res4[metric][linkage]['SoS'].append(SoS)

fi = glob.glob('res4/*')
p = lambda x: x[x.rindex('/')+1:x.rindex('.')].split('_')[1:]
res5 = {}

for f in fi:
    metric, linkage = p(f)
    if not metric in res5:
        res5[metric] = {}
    if not linkage in res5:
        res5[linkage] = {}
    if not linkage in res5[metric]:
        res5[metric][linkage] = {'VI':[],'SoS':[]}
    if not metric in res5[linkage]:
        res5[linkage][metric] = {'VI':[],'SoS':[]}
    open_file = open(f)
    for line in open_file:
        VI, SoS = line.rstrip().split()
        VI = float(VI)
        SoS = float(SoS)
        res5[linkage][metric]['VI'].append(VI)
        res5[linkage][metric]['SoS'].append(SoS)
        res5[metric][linkage]['VI'].append(VI)
        res5[metric][linkage]['SoS'].append(SoS)
        
def get_data(dictionary,metric,linkages,score):
    data = []
    for linkage in linkages:
        to_append = np.array(dictionary[metric][linkage][score])
        to_append.shape = (-1,1)
        data.append(to_append)
    return data

res4_euc_SoS = get_data(res4,'euc',['single','complete','ward'],'SoS')
res5_euc_SoS = get_data(res5,'euc',['single','complete','ward'],'SoS')
res4_rf_SoS = get_data(res4,'rf',['single','complete','ward'],'SoS')
res5_rf_SoS = get_data(res5,'rf',['single','complete','ward'],'SoS')
res4_sym_SoS = get_data(res4,'sym',['single','complete','ward'],'SoS')
res5_sym_SoS = get_data(res5,'sym',['single','complete','ward'],'SoS')
res4_geo_SoS = get_data(res4,'geo',['single','complete','ward'],'SoS')
res5_geo_SoS = get_data(res5,'geo',['single','complete','ward'],'SoS')

res4_euc_VI = get_data(res4,'euc',['single','complete','ward'],'VI')
res5_euc_VI = get_data(res5,'euc',['single','complete','ward'],'VI')
res4_rf_VI = get_data(res4,'rf',['single','complete','ward'],'VI')
res5_rf_VI = get_data(res5,'rf',['single','complete','ward'],'VI')
res4_sym_VI = get_data(res4,'sym',['single','complete','ward'],'VI')
res5_sym_VI = get_data(res5,'sym',['single','complete','ward'],'VI')
res4_geo_VI = get_data(res4,'geo',['single','complete','ward'],'VI')
res5_geo_VI = get_data(res5,'geo',['single','complete','ward'],'VI')

res4_single_SoS = get_data(res4,'single',['sym','rf','euc','geo'],'SoS')
res5_single_SoS = get_data(res5,'single',['sym','rf','euc','geo'],'SoS')
res4_complete_SoS = get_data(res4,'complete',['sym','rf','euc','geo'],'SoS')
res5_complete_SoS = get_data(res5,'complete',['sym','rf','euc','geo'],'SoS')
res4_ward_SoS = get_data(res4,'ward',['sym','rf','euc','geo'],'SoS')
res5_ward_SoS = get_data(res5,'ward',['sym','rf','euc','geo'],'SoS')

res4_single_VI = get_data(res4,'single',['sym','rf','euc','geo'],'VI')
res5_single_VI = get_data(res5,'single',['sym','rf','euc','geo'],'VI')
res4_complete_VI = get_data(res4,'complete',['sym','rf','euc','geo'],'VI')
res5_complete_VI = get_data(res5,'complete',['sym','rf','euc','geo'],'VI')
res4_ward_VI = get_data(res4,'ward',['sym','rf','euc','geo'],'VI')
res5_ward_VI = get_data(res5,'ward',['sym','rf','euc','geo'],'VI')

res4_metrics_VI = res4_euc_VI+res4_rf_VI+res4_sym_VI+res4_geo_VI
res4_metrics_SoS = res4_euc_SoS+res4_rf_SoS+res4_sym_SoS+res4_geo_SoS
res5_metrics_VI = res5_euc_VI+res5_rf_VI+res5_sym_VI+res5_geo_VI
res5_metrics_SoS = res5_euc_SoS+res5_rf_SoS+res5_sym_SoS+res5_geo_SoS

res4_linkages_VI = res4_single_VI+res4_complete_VI+res4_ward_VI
res4_linkages_SoS = res4_single_SoS+res4_complete_SoS+res4_ward_SoS
res5_linkages_VI = res5_single_VI+res5_complete_VI+res5_ward_VI
res5_linkages_SoS = res5_single_SoS+res5_complete_SoS+res5_ward_SoS

fig = figure(figsize=(16,12))
bp = plt.boxplot(res5_linkages_SoS, notch=1, sym='+', vert=1, whis=1.5)
plt.setp(bp['boxes'], color='black')
plt.setp(bp['whiskers'], color='black')
plt.setp(bp['fliers'], color='red', marker='+')
ax1 = fig.add_subplot(111)
ax1.yaxis.grid(True, linestyle='-', color='grey',
              alpha=0.5)
ax1.set_axisbelow(True)
ax1.set_title('Comparison of Performance of Linkage Methods and Distance Metrics',fontsize=16)
ax1.set_xlabel('Distance Metric',fontsize=14)
ax1.set_ylabel('Sum of Squares',fontsize=14)

boxColors = ['lightblue','lightgreen','#FFFFCC']
numBoxes = 12
medians = range(numBoxes)
for i in range(numBoxes):
  box = bp['boxes'][i]
  boxX = []
  boxY = []
  for j in range(10):
      boxX.append(box.get_xdata()[j])
      boxY.append(box.get_ydata()[j])
  boxCoords = zip(boxX,boxY)
  if i in [0,1,2,3]:
      if i==0:
          boxPolygon = Polygon(boxCoords, facecolor=boxColors[0], label="Single")
      else:
          boxPolygon = Polygon(boxCoords, facecolor=boxColors[0])
  elif i in [4,5,6,7]:
      if i == 4:
          boxPolygon = Polygon(boxCoords, facecolor=boxColors[1], label="Complete")
      else:
          boxPolygon = Polygon(boxCoords, facecolor=boxColors[1])
  elif i in [8,9,10,11]:
      if i == 8:
          boxPolygon = Polygon(boxCoords, facecolor=boxColors[2], label="Ward's Method")
      else:
          boxPolygon = Polygon(boxCoords, facecolor=boxColors[2])
  ax1.add_patch(boxPolygon)
  # Now draw the median lines back over what we just filled in
  med = bp['medians'][i]
  medianX = []
  medianY = []
  for j in range(2):
      medianX.append(med.get_xdata()[j])
      medianY.append(med.get_ydata()[j])
      plot(medianX, medianY, 'k')
      medians[i] = medianY[0]
  # Finally, overplot the sample averages, with horixzontal alignment
  # in the center of each box
  plot([np.average(med.get_xdata())], [np.average(res5_linkages_SoS[i])],
           color='r', marker='*', markeredgecolor='k')


legend(title="Linkage Method")
leg = gca().get_legend()
setp(leg.get_texts(), fontsize=12)
xtickNames = setp(ax1, xticklabels=['RF','Weighted RF', 'Felsenstein','Geodesic']*3)
setp(xtickNames, rotation=45, fontsize=12)
fig.savefig('comparison.pdf')
show()

