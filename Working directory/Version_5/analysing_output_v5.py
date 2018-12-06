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

inputFile = 'test_output/output2.tsv'
sampleStats = 'test_output/sample_stats_complete.xlsx'
statsColumn = 'Species'
data = defaultdict(list)
virus = defaultdict(list)
virusLabel = set()
sampleNames = set()
Abs = []
Covs = []
Errs = []
allStats = pd.read_excel(sampleStats)
stats = allStats[['FileName',statsColumn]]

with open(inputFile) as f:
	firstLine = f.readline()  # Don't include header
	for line in f:
		(fileName, virusId, virusName, abundance, coverage, errorRate, percentRemoved) = line.strip().split("\t")
		# if (virusId != 'NC_001422.1') & (virusId != 'NC_022518.1') & (virusId != 'NC_009334.1') & (virusId != 'NC_007605.1') & (virusId != 'NC_001604.1'):
		data[fileName].append((virusId, abundance, coverage, errorRate, percentRemoved))
		virus[virusId].append((fileName, abundance, coverage, errorRate, percentRemoved))
		virusLabel.add(virusId)
		Abs.append(abundance)
		Covs.append(coverage)
		Errs.append(errorRate)
		sampleNames.add(fileName)

outputAbe = defaultdict(list)
outputCov = defaultdict(list)
outputErr = defaultdict(list)
outputPercentRemoved = defaultdict(list)
outputScore = defaultdict(list)
names = []

for sample in data:
	# if (sample != 'sample_data/Pan_troglodytes_troglodytes-13656_Brigitte.mpile.gz'):
		for key,value in virus.items():
				listId = [i[0] for i in value]
				listAbe = [j[1] for j in value]
				listCov = [j[2] for j in value]
				listErr = [j[3] for j in value]
				listPerRem = [j[4] for j in value]
				if any(elem == sample for elem in listId):
					ind = listId.index(sample)
					outputAbe[sample].append(listAbe[ind])
					outputCov[sample].append(listCov[ind])
					outputErr[sample].append(listErr[ind])
					outputPercentRemoved[sample].append(listPerRem[ind])
					outputScore[sample].append(round(float(listAbe[ind])*float(listCov[ind])*(1-float(listErr[ind]))*(1-float(listPerRem[ind])),4))
				else:
					outputAbe[sample].append(0)  # If the individual doesn't have the virus add a 0 instead
					outputCov[sample].append(0)
					outputErr[sample].append(0)
					outputScore[sample].append(0)

# for name in sampleNames:
# 	print stats.loc[stats.itemsets == (name,), 'Source']

# SCORE
# Create pandas DataFrame for Super score
outScore = pd.DataFrame.from_dict(outputScore, orient='index', columns=virusLabel)  # Convert the dictionary to a DataFrame
index = outScore.index  # Sample rows
cols = list(virusLabel)  # Virus columns
valScore = outScore.values
valScore_NaN = np.ma.masked_where(valScore == 0, valScore)  # Abundance == 0 is replaced by NaN
dfScore = pd.DataFrame(valScore_NaN, index=index, columns=cols)
dfScore = dfScore[dfScore.columns].astype(float).sort_index()
# print dfScore
#
# dfScoreBinary = pd.DataFrame(valScore, index=index, columns=cols)
# dfScoreBinary = dfScoreBinary[dfScoreBinary.columns].astype(float).sort_index()
# Zero = dfScoreBinary.where(dfScoreBinary != 0)
# print Zero
# for ids in virusLabel:
# 	for index, row in dfScoreBinary.iterrows():
# 		print dfScoreBinary[row[ids]]
# 		if row[ids] < 5:
# 			dfScoreBinary[index] = 0
# 		else:
# 			dfScoreBinary[index] = 1
# print dfScoreBinary

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
# # ERROR RATE
# # Create pandas DataFrame for Coverage
# outErr = pd.DataFrame.from_dict(outputErr, orient='index', columns=virusLabel)  # Convert the dictionary to a DataFrame
# valErr = outErr.values
# valErr_NaN = np.ma.masked_where(valErr == 0, valErr)  # Error rate == 0 is replaced by NaN
# dfErr = pd.DataFrame(valErr_NaN, index=index, columns=cols)
# dfErr = dfErr[dfErr.columns].astype(float)

# The DataFrames will have the format:
#                                                    Rose rosette emaravirus  ...   Bwamba orthobunyavirus
# ../../sample_data/Pan_troglodytes_verus-N016_Al...            NaN							NaN
# ../../sample_data/Homo_sapiens_C07.mpile.gz                   16.22 						NaN
# ../../sample_data/Pan_troglodytes_schweinfurthi...            NaN							NaN
#  ...
# ../../sample_data/Pongo_abelii-A948_Kiki.mpile.gz 			NaN							6.37
# ../../sample_data/Pongo_abelii-A948.mpile.gz 					NaN							4.34

# HEAT MAPS

# General color scheme for the heat maps
cmap = plt.cm.Greens
cmap.set_bad(color='white')
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

# SUPER SCORE
plt.figure(1, figsize = (8,8))
plt.yticks(np.arange(len(dfScore.index)), dfScore.index)  # Set sample names as Y-axis
plt.xticks(np.arange(len(dfScore.columns)), dfScore.columns, rotation='vertical')  # Set virus ids as X-axis
ax = plt.gca()  # Get the axes
heatMapScore = plt.imshow(dfScore, cmap=cmap)
divider = make_axes_locatable(ax)  # Adjust the dimensions of the color bar
cax = divider.append_axes("right", size="3%", pad=0.05)
plt.suptitle('Super Score')
plt.colorbar(heatMapScore, cax)

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
# COVERAGE
# Heat map using output values for Coverage
# plt.figure(3)
# plt.yticks(np.arange(len(dfCov.index)), dfCov.index)  # Set sample names as Y-axis
# plt.xticks(np.arange(len(dfCov.columns)), dfCov.columns, rotation='vertical')  # Set virus ids as X-axis
# ax = plt.gca()  # Get the axes
# heatMapCov = plt.imshow(dfCov, cmap=cmap)
# divider = make_axes_locatable(ax)  # Adjust the dimensions of the color bar
# cax = divider.append_axes("right", size="3%", pad=0.05)
# plt.suptitle('Coverage')
# plt.colorbar(heatMapCov, cax)
# plt.show()
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
#
# vars = np.column_stack((np.asarray(Abs).astype(np.float), np.asarray(Covs).astype(np.float), np.asarray(Errs).astype(np.float)))
# print np.corrcoef(vars, rowvar=False)
#
# fig = plt.figure(7)
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(np.asarray(Abs).astype(np.float), np.asarray(Covs).astype(np.float), np.asarray(Errs).astype(np.float), c='r')
# ax.set_xlabel('Abundance')
# ax.set_ylabel('Coverage')
# ax.set_zlabel('Error Rate')
#
# PCA

outPCA = pd.DataFrame.from_dict(outputScore, orient='index')  # Convert the dictionary to a DataFrame
index = outPCA.index  # Sample rows
valPCA = outPCA.values
dfPCA = pd.DataFrame(valPCA, index=index)
dfPCA = dfPCA[dfPCA.columns].astype(float)

statsDict = defaultdict(list)
for index in dfPCA.index.tolist():
	for name in stats['FileName']:
		if index != 'Pan_troglodytes_troglodytes-13656_Brigitte':
			if index == name:
				statsDict[name].append(stats.at[stats[stats['FileName'] == name].index.values.astype(int)[0], statsColumn])
# print statsDict
dfStats = pd.DataFrame.from_dict(statsDict, orient='index', columns=[statsColumn])  # Convert the dictionary to a DataFrame
X = StandardScaler().fit_transform(dfPCA)

pca = PCA(n_components=2)
principalComponents = pca.fit_transform(X)
principalDf = pd.DataFrame(data=principalComponents, columns=['principal component 1', 'principal component 2'], index=sampleNames)
# print principalDf['principal component 2']
# dfs = [principalDf, dfStats]
# colsPCA = ['principal component 1', 'principal component 2', statsColumn]
# keys = ['value1', 'value2']
# pd.concat(
#     [df.set_index(colsPCA).value for df in dfs],
#     axis=1, keys=keys)
# print dfs

finalDf = pd.concat([principalDf, dfStats], axis=1, sort=False)
# print finalDf['principal component 2']
print principalDf.loc[principalDf['principal component 2'] == max(principalDf['principal component 2'])]
# print principalDf.idxmax()
# print principalDf['principal component 2'][76]
fig = plt.figure(2,figsize = (8,8))
ax = fig.add_subplot(1,1,1)
ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)
ax.set_title('2 component PCA', fontsize= 20)
targets = ['Gorilla', 'Homo sapiens', 'Pan troglodytes', 'Pan paniscus', 'Pongo']
colors = ['r','g','b','m','y']
for target, color in zip(targets,colors):
    indicesToKeep = finalDf[statsColumn] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 2']
               , c = color
               , s = 50)
ax.legend(targets)
ax.grid()
plt.show()