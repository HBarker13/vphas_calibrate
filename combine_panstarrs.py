#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python


""" Combine the four panstarrs files needed to cover a vphas pointing. Remove overlaps """

import os
import csv
import glob
import itertools

#panstarrs files
pan_names = [l for l in glob.glob( os.getcwd()+'/panstarrs*.csv') if 'combined' not in l]


combined_list = []
for name in pan_names:
	print 'Reading in ', name
	with open(name, 'rb') as f:
		reader = csv.reader(f)
		for ind,line in enumerate(reader):
			if ind==0: continue #skip column names
			combined_list.append(line)

print
combined_list.sort()
combined = list( k for k,_ in itertools.groupby(combined_list) )


print combined[0]
from matplotlib import pyplot as plt

x = [ float(l[2]) for l in combined ] 
print x[0:5]
y = [ float(l[3]) for l in combined ]


plt.plot(x,y, 'o')
plt.show()




