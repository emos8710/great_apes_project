from collections import defaultdict
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA as PCA
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

"""
Script that produces heat maps from the output file from main.py. 
"""

inputFile = 'test_output/output.tsv'
data = defaultdict(list)
virus = defaultdict(list)
virusLabel = set()

Abs = []
Covs = []
Errs = []

with open(inputFile) as f:
	firstLine = f.readline()  # Don't include header
	for line in f:
		(fileName, virusId, virusName, abundance, coverage, errorRate) = line.strip().split("\t")
		if (virusId != 'NC_001422.1') & (virusId != 'NC_022518.1') & (virusId != 'NC_009334.1') & (virusId != 'NC_007605.1') & (virusId != 'NC_001604.1'):
			data[fileName].append((virusId, abundance, coverage, errorRate))
			virus[virusId].append((fileName, abundance, coverage, errorRate))
			virusLabel.add(virusId)
			Abs.append(abundance)
			Covs.append(coverage)
			Errs.append(errorRate)

outputAbe = defaultdict(list)
outputCov = defaultdict(list)
outputErr = defaultdict(list)
apes = {}
for sample in data:
	score = []
	for key, value in virus.items():
		listId = [i[0] for i in value]
		listAbe = [j[1] for j in value]
		listCov = [j[2] for j in value]
		listErr = [j[3] for j in value]
		if any(elem == sample for elem in listId):
			ind = listId.index(sample)
			outputAbe[sample].append(listAbe[ind])
			outputCov[sample].append(listCov[ind])
			outputErr[sample].append(listErr[ind])
		else:
			outputAbe[sample].append(0)  # If the individual doesn't have the virus add a 0 instead
			outputCov[sample].append(0)
			outputErr[sample].append(0)
 		score.append((float(outputAbe[sample][-1])/100) + 2*float(outputCov[sample][-1]) + 2*float(outputErr[sample][-1]))
	apes[sample] = score
df = pd.DataFrame.from_dict(apes, orient='index', columns=virusLabel)  # Convert the dictionary to a DataFrame
x = StandardScaler().fit_transform(df)

pca = PCA(n_components=3)
principalComponents = pca.fit_transform(x)
principalDf = pd.DataFrame(data=principalComponents, columns=['principal component 1', 'principal component 2', 'principal component 3'])

fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection='3d')
ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)
ax.set_title('2 component PCA', fontsize = 20)

xaxis = []
yaxis = []
zaxis = []
for i in range(len(principalDf.values)):
	xaxis.append(principalDf.values[i][0])
	yaxis.append(principalDf.values[i][1])
	zaxis.append(principalDf.values[i][2])
ax.scatter(xaxis, yaxis, zaxis)
plt.show()