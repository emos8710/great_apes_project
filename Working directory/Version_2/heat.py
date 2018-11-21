from collections import defaultdict
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

"""
Preliminary script that produces a heatmap (or more) from the output file from main.py. So far this is only done for
abundance. 
"""

inputFile = '../../test_output/output.tsv'
data = defaultdict(list)
virus = defaultdict(list)
virusLabel = set()
sampleNames = set()

Abs = []
Covs = []
Errs = []

with open(inputFile) as f:
	firstLine = f.readline()  # Don't include header
	for line in f:
		(fileName, virusId, virusName, abundance, coverage, errorRate) = line.strip().split("\t")
		data[fileName].append((virusName, abundance))
		virus[virusName].append((fileName, abundance))
		virusLabel.add(virusName)
		sampleNames.add(fileName)
		Abs.append(abundance)
		Covs.append(coverage)
		Errs.append(errorRate)
f.close()

output = defaultdict(list)  # Store which samples contain which viruses in a dictionary (based on abundance)
for sample in data:
	for key,value in virus.items():
		listNames = [i[0] for i in value]
		listValues = [j[1] for j in value]
		if any(elem == sample for elem in listNames):
			ind = listNames.index(sample)
			output[sample].append(listValues[ind])
		else:
			output[sample].append(0)  # If the individual doesn't have the virus


pdOut = pd.DataFrame.from_dict(output, orient='index', columns=virusLabel)  # Convert the dictionary to a DataFrame

# The DataFrame will have the format:
#                                                    Rose rosette emaravirus  ...   Bwamba orthobunyavirus
# ../../sample_data/Pan_troglodytes_verus-N016_Al...            0							0
# ../../sample_data/Homo_sapiens_C07.mpile.gz                   16.22 						0
# ../../sample_data/Pan_troglodytes_schweinfurthi...            0							0
#  ...
# ../../sample_data/Pongo_abelii-A948_Kiki.mpile.gz 			0							6.37
# ../../sample_data/Pongo_abelii-A948.mpile.gz 					0							4.34

index = pdOut.index
cols = list(virusLabel)
val = pdOut.values

# val = np.ma.masked_where(val == 0, val)  # Set color to white for missing viruses (abundance == 0)
# cmap = plt.cm.Greens
# cmap.set_bad(color='white')
#
# df = pd.DataFrame(val, index=index, columns=cols)  # Make another DataFrame (lol) this just made it easier to plot.
# df = df[df.columns].astype(float)  # Change dtype to float64
# plt.yticks(np.arange(len(df.index)), df.index)  # Set sample names as Y-axis
# plt.xticks(np.arange(len(df.columns)), df.columns, rotation='vertical')  # Set virus names as X-axis
#
# ax = plt.gca()  # Get the axes
# heatMap = plt.imshow(df, cmap=cmap)
# divider = make_axes_locatable(ax)  # Adjust the dimensions of the colorbar
# cax = divider.append_axes("right", size="3%", pad=0.05)
# plt.colorbar(heatMap, cax)
# plt.show()

# CORRELATION

vars = np.column_stack((np.asarray(Abs).astype(np.float), np.asarray(Covs).astype(np.float), np.asarray(Errs).astype(np.float)))
print vars
print np.corrcoef(vars, rowvar=False)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(np.asarray(Abs).astype(np.float), np.asarray(Covs).astype(np.float), np.asarray(Errs).astype(np.float), c='r')
ax.set_xlabel('Abundance')
ax.set_ylabel('Coverage')
ax.set_zlabel('Error Rate')
plt.show()
