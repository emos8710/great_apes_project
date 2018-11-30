from collections import defaultdict
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

"""
Script that produces heat maps, 3D scatter-plots and PCA using the output file from main.py. 
"""

# Function which calculates Z-score
def zscore(dataFrame):
	dataFrame = dataFrame[dataFrame.columns].astype(float)
	for col in list(dataFrame.columns):
		dataFrame[col] = (dataFrame[col] - dataFrame[col].mean())/dataFrame[col].std(ddof=1)
	return dataFrame

inputFile = '../../test_output/output.tsv'
sampleStats = '../../test_output/sample_stats_complete.xlsx'
data = defaultdict(list)
virus = defaultdict(list)
virusLabel = set()
sampleNames = []
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
			sampleNames.append(fileName)

outputAbe = defaultdict(list)
outputCov = defaultdict(list)
outputErr = defaultdict(list)
outputPCA = defaultdict(list)
names = []

for sample in data:
	if (sample != 'sample_data/Pan_troglodytes_schweinfurthii-Mgbadolite.mpile.gz') & (sample != 'sample_data/Gorilla_beringei_graueri-A929_Kaisi.mpile.gz') & ( sample != 'sample_data/Pongo_pygmaeus-A946.mpile.gz'):
		for key,value in virus.items():
				listId = [i[0] for i in value]
				listAbe = [j[1] for j in value]
				listCov = [j[2] for j in value]
				listErr = [j[3] for j in value]
				if any(elem == sample for elem in listId):
					ind = listId.index(sample)
					outputAbe[sample].append(listAbe[ind])
					outputCov[sample].append(listCov[ind])
					outputErr[sample].append(listErr[ind])
					outputPCA[sample].append(listAbe[ind])
					outputPCA[sample].append(listCov[ind])
					outputPCA[sample].append(listErr[ind])
				else:
					outputAbe[sample].append(0)  # If the individual doesn't have the virus add a 0 instead
					outputCov[sample].append(0)
					outputErr[sample].append(0)
					outputPCA[sample].append(0)
					outputPCA[sample].append(0)
					outputPCA[sample].append(0)

# ABUNDANCE
# Create pandas DataFrame for Abundance
outAbe = pd.DataFrame.from_dict(outputAbe, orient='index', columns=virusLabel)  # Convert the dictionary to a DataFrame
index = outAbe.index  # Sample rows
cols = list(virusLabel)  # Virus columns
valAbe = outAbe.values
valAbe_NaN = np.ma.masked_where(valAbe == 0, valAbe)  # Abundance == 0 is replaced by NaN
dfAbe = pd.DataFrame(valAbe_NaN, index=index, columns=cols)
dfAbe = dfAbe[dfAbe.columns].astype(float)
# COVERAGE
# Create pandas DataFrame for Coverage
outCov = pd.DataFrame.from_dict(outputCov, orient='index', columns=virusLabel)  # Convert the dictionary to a DataFrame
valCov = outCov.values
valCov_NaN = np.ma.masked_where(valCov == 0, valCov)  # Coverage == 0 is replaced by NaN
dfCov = pd.DataFrame(valCov_NaN, index=index, columns=cols)
dfCov = dfCov[dfCov.columns].astype(float)
# ERROR RATE
# Create pandas DataFrame for Coverage
outErr = pd.DataFrame.from_dict(outputErr, orient='index', columns=virusLabel)  # Convert the dictionary to a DataFrame
valErr = outErr.values
valErr_NaN = np.ma.masked_where(valErr == 0, valErr)  # Error rate == 0 is replaced by NaN
dfErr = pd.DataFrame(valErr_NaN, index=index, columns=cols)
dfErr = dfErr[dfErr.columns].astype(float)

# The DataFrames will have the format:
#                                                    Rose rosette emaravirus  ...   Bwamba orthobunyavirus
# ../../sample_data/Pan_troglodytes_verus-N016_Al...            NaN							NaN
# ../../sample_data/Homo_sapiens_C07.mpile.gz                   16.22 						NaN
# ../../sample_data/Pan_troglodytes_schweinfurthi...            NaN							NaN
#  ...
# ../../sample_data/Pongo_abelii-A948_Kiki.mpile.gz 			NaN							6.37
# ../../sample_data/Pongo_abelii-A948.mpile.gz 					NaN							4.34

# HEAT MAPS
# #
# # # General color scheme for the heat maps
# cmap = plt.cm.Greens
# cmap.set_bad(color='white')
#
# # ABUNDANCE
# # Heat map using output values for Abundance
# plt.figure(1)
# plt.yticks(np.arange(len(dfAbe.index)), dfAbe.index)  # Set sample names as Y-axis
# plt.xticks(np.arange(len(dfAbe.columns)), dfAbe.columns, rotation='vertical')  # Set virus ids as X-axis
# ax = plt.gca()  # Get the axes
# heatMapAbe = plt.imshow(dfAbe, cmap=cmap)
# divider = make_axes_locatable(ax)  # Adjust the dimensions of the color bar
# cax = divider.append_axes("right", size="3%", pad=0.05)
# plt.suptitle('Abundance')
# plt.colorbar(heatMapAbe, cax)
#
# # Heat map using z-score for Abundance
# dfAbe0 = pd.DataFrame(valAbe, index=index, columns=cols)
# dfAbeZ = zscore(dfAbe0)
# plt.figure(2)
# plt.yticks(np.arange(len(dfAbeZ.index)), dfAbeZ.index)  # Set sample names as Y-axis
# plt.xticks(np.arange(len(dfAbeZ.columns)), dfAbeZ.columns, rotation='vertical')  # Set virus ids as X-axis
# ax = plt.gca()  # Get the axes
# heatMapAbeZ = plt.imshow(dfAbeZ, cmap=cmap)
# divider = make_axes_locatable(ax)  # Adjust the dimensions of the color bar
# cax = divider.append_axes("right", size="3%", pad=0.05)
# plt.suptitle('Abundance Z-score')
# plt.colorbar(heatMapAbeZ, cax)
#
# # COVERAGE
# # Heat map using output values for Coverage
# plt.figure(3)
# plt.yticks(np.arange(len(dfCov.index)), dfCov.index)  # Set sample names as Y-axis
# plt.xticks(np.arange(len(dfCov.columns)), dfCov.columns, rotation='vertical')  # Set virus ids as X-axis
# ax = plt.gca()  # Get the axes
# heatMapCov = plt.imshow(dfCov, cmap=cmap)
# divider = make_axes_locatable(ax)  # Adjust the dimensions of the color bar
# cax = divider.append_axes("right", size="3%", pad=0.05)
# plt.suptitle('Coverage')
# plt.colorbar(heatMapCov, cax)
#
# # Heat map using z-score for Coverage
# dfCov0 = pd.DataFrame(valCov, index=index, columns=cols)
# dfCovZ = zscore(dfCov0)
# plt.figure(4)
# plt.yticks(np.arange(len(dfCovZ.index)), dfCovZ.index)  # Set sample names as Y-axis
# plt.xticks(np.arange(len(dfCovZ.columns)), dfCovZ.columns, rotation='vertical')  # Set virus ids as X-axis
# ax = plt.gca()  # Get the axes
# heatMapCovZ = plt.imshow(dfCovZ, cmap=cmap)
# divider = make_axes_locatable(ax)  # Adjust the dimensions of the color bar
# cax = divider.append_axes("right", size="3%", pad=0.05)
# plt.suptitle('Coverage Z-score')
# plt.colorbar(heatMapCovZ, cax)
#
# # ERROR RATE
# # Heat map using output values for Error rate
# plt.figure(5)
# plt.yticks(np.arange(len(dfErr.index)), dfErr.index)  # Set sample names as Y-axis
# plt.xticks(np.arange(len(dfErr.columns)), dfErr.columns, rotation='vertical')  # Set virus ids as X-axis
# ax = plt.gca()  # Get the axes
# heatMapErr = plt.imshow(dfErr, cmap=cmap)
# divider = make_axes_locatable(ax)  # Adjust the dimensions of the color bar
# cax = divider.append_axes("right", size="3%", pad=0.05)
# plt.suptitle('Error rate')
# plt.colorbar(heatMapErr, cax)
#
# # Heat map using z-score for Error rate
# dfErr0 = pd.DataFrame(valErr, index=index, columns=cols)
# dfErrZ = zscore(dfErr0)
# plt.figure(6)
# plt.yticks(np.arange(len(dfErrZ.index)), dfErrZ.index)  # Set sample names as Y-axis
# plt.xticks(np.arange(len(dfErrZ.columns)), dfErrZ.columns, rotation='vertical')  # Set virus ids as X-axis
# ax = plt.gca()  # Get the axes
# heatMapErrZ = plt.imshow(dfErrZ, cmap=cmap)
# divider = make_axes_locatable(ax)  # Adjust the dimensions of the color bar
# cax = divider.append_axes("right", size="3%", pad=0.05)
# plt.suptitle('Error rate Z-score')
# plt.colorbar(heatMapErrZ, cax)

# CORRELATION

vars = np.column_stack((np.asarray(Abs).astype(np.float), np.asarray(Covs).astype(np.float), np.asarray(Errs).astype(np.float)))
print np.corrcoef(vars, rowvar=False)

fig = plt.figure(7)
ax = fig.add_subplot(111, projection='3d')
ax.scatter(np.asarray(Abs).astype(np.float), np.asarray(Covs).astype(np.float), np.asarray(Errs).astype(np.float), c='r')
ax.set_xlabel('Abundance')
ax.set_ylabel('Coverage')
ax.set_zlabel('Error Rate')

# PCA
stats = pd.read_excel(sampleStats)

#
# virusList = []
# for key in virus:
# 	virusList.append((virus[key][0][1],virus[key][0][2],virus[key][0][3],virus[key][0][0]))
#
# outPCA = pd.DataFrame.from_dict(outputPCA, orient='index')  # Convert the dictionary to a DataFrame
# index = outPCA.index  # Sample rows
# valPCA = outPCA.values
# dfPCA = pd.DataFrame(valPCA, index=index)
# dfPCA = dfPCA[dfPCA.columns].astype(float)
# # print dfPCA.iloc[[93]]
# X = StandardScaler().fit_transform(dfPCA)
#
# pca = PCA(n_components=2)
# principalComponents = pca.fit_transform(X)
# principalDf = pd.DataFrame(data=principalComponents, columns=['principal component 1', 'principal component 2'])
# # print principalDf.loc[principalDf['principal component 2'] == max(principalDf['principal component 2'])]
# # print principalDf.idxmax()
# # print principalDf['principal component 2'][76]
# fig = plt.figure(8,figsize = (8,8))
# ax = fig.add_subplot(1,1,1)
# ax.set_xlabel('Principal Component 1', fontsize = 15)
# ax.set_ylabel('Principal Component 2', fontsize = 15)
# ax.set_title('2 component PCA', fontsize= 20)
# ax.scatter(principalDf['principal component 1'], principalDf['principal component 2'])
# ax.grid()
# plt.show()
