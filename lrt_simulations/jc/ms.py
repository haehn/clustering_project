import glob

files = glob.glob('./clustersuccess/fullresults/*.txt')
l=[]
for i in [2,4,6,8,10,12]:
    for j in range(1,101):
        l.append('./clustersuccess/fullresults/rf{0}_{1}.txt'.format(i,j))

for file in files:
    l.remove(file)

print l

