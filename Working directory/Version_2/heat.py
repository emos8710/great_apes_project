from collections import defaultdict
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA as sklearnPCA
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

"""
Script that produces heat maps from the output file from main.py. 
"""

# Function which calculates Z-score
def zscore(dataFrame):
	dataFrame = dataFrame[dataFrame.columns].astype(float)
	for col in list(dataFrame.columns):
		dataFrame[col] = (dataFrame[col] - dataFrame[col].mean())/dataFrame[col].std(ddof=1)
	return dataFrame

inputFile = '../../test_output/output.tsv'
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

for sample in data:
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
			else:
				outputAbe[sample].append(0)  # If the individual doesn't have the virus add a 0 instead
				outputCov[sample].append(0)
				outputErr[sample].append(0)

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
#
# # General color scheme for the heat maps
# cmap = plt.cm.Greens
# cmap.set_bad(color='white')
#
# # ABUNDANCE
# # Heat map using output values for Abundance
# plt.yticks(np.arange(len(dfAbe.index)), dfAbe.index)  # Set sample names as Y-axis
# plt.xticks(np.arange(len(dfAbe.columns)), dfAbe.columns, rotation='vertical')  # Set virus ids as X-axis
# ax = plt.gca()  # Get the axes
# heatMapAbe = plt.imshow(dfAbe, cmap=cmap)
# divider = make_axes_locatable(ax)  # Adjust the dimensions of the color bar
# cax = divider.append_axes("right", size="3%", pad=0.05)
# plt.suptitle('Abundance')
# plt.colorbar(heatMapAbe, cax)
# plt.show()
# # Heat map using z-score for Abundance
# dfAbe0 = pd.DataFrame(valAbe, index=index, columns=cols)
# dfAbeZ = zscore(dfAbe0)
# plt.yticks(np.arange(len(dfAbeZ.index)), dfAbeZ.index)  # Set sample names as Y-axis
# plt.xticks(np.arange(len(dfAbeZ.columns)), dfAbeZ.columns, rotation='vertical')  # Set virus ids as X-axis
# heatMapAbeZ = plt.imshow(dfAbeZ, cmap=cmap)
# plt.suptitle('Abundance Z-score')
# plt.colorbar(heatMapAbeZ, cax)
# plt.show()
#
# # COVERAGE
# # Heat map using output values for Coverage
# plt.yticks(np.arange(len(dfCov.index)), dfCov.index)  # Set sample names as Y-axis
# plt.xticks(np.arange(len(dfCov.columns)), dfCov.columns, rotation='vertical')  # Set virus ids as X-axis
# heatMapCov = plt.imshow(dfCov, cmap=cmap)
# plt.suptitle('Coverage')
# plt.colorbar(heatMapCov, cax)
# plt.show()
#
# # Heat map using z-score for Coverage
# dfCov0 = pd.DataFrame(valCov, index=index, columns=cols)
# dfCovZ = zscore(dfCov0)
# plt.yticks(np.arange(len(dfCovZ.index)), dfCovZ.index)  # Set sample names as Y-axis
# plt.xticks(np.arange(len(dfCovZ.columns)), dfCovZ.columns, rotation='vertical')  # Set virus ids as X-axis
# heatMapCovZ = plt.imshow(dfCovZ, cmap=cmap)
# plt.suptitle('Coverage Z-score')
# plt.colorbar(heatMapCovZ, cax)
# plt.show()
#
# # ERROR RATE
# # Heat map using output values for Error rate
# plt.yticks(np.arange(len(dfErr.index)), dfErr.index)  # Set sample names as Y-axis
# plt.xticks(np.arange(len(dfErr.columns)), dfErr.columns, rotation='vertical')  # Set virus ids as X-axis
# heatMapErr = plt.imshow(dfErr, cmap=cmap)
# plt.suptitle('Error rate')
# plt.colorbar(heatMapErr, cax)
# plt.show()
#
# # Heat map using z-score for Error rate
# dfErr0 = pd.DataFrame(valErr, index=index, columns=cols)
# dfErrZ = zscore(dfErr0)
# plt.yticks(np.arange(len(dfErrZ.index)), dfErrZ.index)  # Set sample names as Y-axis
# plt.xticks(np.arange(len(dfErrZ.columns)), dfErrZ.columns, rotation='vertical')  # Set virus ids as X-axis
# heatMapErrZ = plt.imshow(dfErrZ, cmap=cmap)
# plt.suptitle('Error rate Z-score')
# plt.colorbar(heatMapErrZ, cax)
# plt.show()
#
# # CORRELATION
#
# vars = np.column_stack((np.asarray(Abs).astype(np.float), np.asarray(Covs).astype(np.float), np.asarray(Errs).astype(np.float)))
# print np.corrcoef(vars, rowvar=False)
#
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(np.asarray(Abs).astype(np.float), np.asarray(Covs).astype(np.float), np.asarray(Errs).astype(np.float), c='r')
# ax.set_xlabel('Abundance')
# ax.set_ylabel('Coverage')
# ax.set_zlabel('Error Rate')
# plt.show()

# for key in virus:
# 	print key, virus[key][0][0], virus[key][0][1], virus[key][0][2], virus[key][0][3]

# PCA

virusList = []
for key in virus:
	virusList.append((virus[key][0][1],virus[key][0][2],virus[key][0][3],virus[key][0][0]))

