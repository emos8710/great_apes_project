from collections import defaultdict
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

"""
Script that produces heat maps, 3D scatter plot and PCA using the output file from main.py and a excel file 
which contains sample statistics (sample_stats.xlsx). The visualization is done by combining the four output columns, 
Abundance, Mapped Percent, Error rate and Procent Removed,
by multiplying them with each other. 
"""

inputFile = 'test_output/output.tsv'
data = defaultdict(list)
virus = defaultdict(list)
virusLabel = set()
sampleNames = set()

sampleStats = 'test_output/sample_stats_complete.xlsx'
statsColumn = 'Birth origin' # Specify which column should be used for PCA
# For ANOVA
# allStats = pd.read_excel(sampleStats, index_col='FileName')
# stats = allStats[statsColumn]

allStats = pd.read_excel(sampleStats)
stats = allStats[['FileName', statsColumn]]
# print stats
Abn = []  # Store for correlation test and 3D scatter plot
Cov = []
Err = []

# Highly present and uninteresting virus to filter away
virus1 = 'NC_001422.1'  # phiX174
virus2 = 'NC_022518.1'  # Human endogenous retrovirus K113
virus3 = 'NC_009334.1'  # Human herpesvirus 4
virus4 = 'NC_007605.1'  # Human gammaherpesvirus 4
virus5 = 'NC_001604.1'  # Enterobacteria phage T7

# If you want to remove outliers from PCA
outlier1 = 'Pongo_pygmaeus-A988'
outlier2 = 'Gorilla_gorilla_gorilla-KB5792_Carolyn'
outlier3 = 'Gorilla_gorilla_gorilla-KB7973_Porta'
with open(inputFile) as f:
	firstLine = f.readline()  # Don't include header
	for line in f:
		(fileName, virusId, virusName, abundance, coverage, errorRate, percentRemoved) = line.strip().split("\t")
		# if (virusId != virus1) & (virusId != virus2) & (virusId != virus3) & (virusId != virus4) & (virusId != virus5):
		# if (fileName != outlier1) & (fileName != outlier2) & (fileName != outlier3):
		data[fileName].append((virusId, abundance, coverage, errorRate, percentRemoved))
		virus[virusId].append((fileName, abundance, coverage, errorRate, percentRemoved))
		print virusId, virusName
		virusLabel.add(virusId)
		sampleNames.add(fileName)
		Abn.append(abundance)
		Cov.append(coverage)
		Err.append(errorRate)

outputScore = defaultdict(list)
names = []

for sample in data:
		for key,value in virus.items():
				listId = [i[0] for i in value]
				listAbe = [j[1] for j in value]
				listCov = [j[2] for j in value]
				listErr = [j[3] for j in value]
				listPerRem = [j[4] for j in value]
				if any(elem == sample for elem in listId):
					ind = listId.index(sample)
					outputScore[sample].append(round(float(listAbe[ind])*float(listCov[ind])* # Multiply values to make a Super score
										(1-float(listErr[ind]))*(1-float(listPerRem[ind])),4))
				else:
					outputScore[sample].append(0)  # If the individual doesn't have the virus a zero is added

# SUPER SCORE
# Create pandas DataFrame for Super score
outScore = pd.DataFrame.from_dict(outputScore, orient='index', columns=virusLabel)  # Convert the dictionary to a DataFrame
index = outScore.index  # Sample rows
cols = list(virusLabel)  # Virus columns
valScore = outScore.values
valScore_NaN = np.ma.masked_where(valScore == 0, valScore)  # Score == 0 is replaced by NaN
dfScore = pd.DataFrame(valScore_NaN, index=index, columns=cols)
dfScore = dfScore[dfScore.columns].astype(float).sort_index()
dfScore0 = pd.DataFrame(valScore, index=index, columns=cols)  # Score == 0 kept as 0
dfScore0 = dfScore0[dfScore0.columns].astype(float).sort_index()
######## NOT DONE

def binary(df):
	scoreThreshold = 0.5  # Set threshold
	for i in df.index:
		for col in df.columns:
			if df.at[i, col] >= scoreThreshold:
				df.at[i,col] = 1
			else:
				df.at[i, col] = 0
	return df

dfScoreBinary = pd.DataFrame(valScore, index=index, columns=cols)  # Score == 0 kept as 0
dfScoreBinary = dfScoreBinary[dfScoreBinary.columns].astype(float).sort_index()
binary(dfScoreBinary)

########

writer = pd.ExcelWriter('test_output/test_data_diff_scores_file.xlsx', engine='xlsxwriter')
dfScore.to_excel(writer, sheet_name='Test_data_scores')
writer.save()

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

# # SUPER SCORE
plt.figure(1, figsize = (8,8))
plt.yticks(np.arange(len(dfScore.index)), dfScore.index)  # Set sample names as Y-axis
plt.xticks(np.arange(len(dfScore.columns)), dfScore.columns, rotation='vertical')  # Set virus ids as X-axis
ax = plt.gca()  # Get the axes
heatMapScore = plt.imshow(dfScore, cmap=cmap)
divider = make_axes_locatable(ax)  # Adjust the dimensions of the color bar
cax = divider.append_axes("right", size="3%", pad=0.05)
plt.suptitle('Super Score')
plt.colorbar(heatMapScore, cax)

# BINARY SCORE
plt.figure(2, figsize = (8,8))
plt.yticks(np.arange(len(dfScoreBinary.index)), dfScoreBinary.index)  # Set sample names as Y-axis
plt.xticks(np.arange(len(dfScoreBinary.columns)), dfScoreBinary.columns, rotation='vertical')  # Set virus ids as X-axis
heatMapScoreBinary = plt.imshow(dfScoreBinary, cmap=cmap)
plt.suptitle('Binary Score')
plt.show()
CORRELATION
vars = np.column_stack((np.asarray(Abn).astype(np.float), np.asarray(Cov).astype(np.float), np.asarray(Err).astype(np.float)))
print np.corrcoef(vars, rowvar=False)

fig = plt.figure(7)
ax = fig.add_subplot(111, projection='3d')
ax.scatter(np.asarray(Abn).astype(np.float), np.asarray(Cov).astype(np.float), np.asarray(Err).astype(np.float), c='r')
ax.set_xlabel('Abundance')
ax.set_ylabel('Coverage')
ax.set_zlabel('Error Rate')

# PCA
statsDict = defaultdict(list)
for index in dfScoreBinary.index.tolist():
	for name in stats['FileName']:
			if index == name:
				# print statsDict[name].append(stats.at[name])
				statsDict[name].append(stats.at[stats[stats['FileName'] == name].index.values.astype(int)[0], statsColumn])

dfStats = pd.DataFrame.from_dict(statsDict, orient='index', columns=[statsColumn])  # Convert the dictionary to a DataFrame
print dfStats
X = StandardScaler().fit_transform(dfScoreBinary)
# print dfStats
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(X)
principalDf = pd.DataFrame(data=principalComponents, columns=['principal component 1', 'principal component 2'], index=sampleNames)
finalDf = pd.concat([principalDf, dfStats], axis=1, sort=False)
print finalDf
# Check for outliers
print principalDf.loc[principalDf['principal component 2'] == max(principalDf['principal component 2'])]

# Plot PCA
fig = plt.figure(3,figsize = (8,8))
ax = fig.add_subplot(1,1,1)
ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)
ax.set_title('2 component PCA', fontsize= 20)
targets = ['Wild born', 'Captive born', 'Human']  # Change according to values for chosen column
colors = ['r','b','g']  # Change according to values for chosen column
for target, color in zip(targets,colors):
	indicesToKeep = finalDf[statsColumn] == target
	ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 2']
               , c = color
               , s = 50)
ax.legend(targets)
ax.grid()
plt.show()

# ANOVA
# ind = []
# for index in dfScore0.index:
# 	for filename in allStats.index:
# 		if index == filename:
# 			ind.append(index)
#
# # print allStats.ix[ind]
# # print dfScore0
# dfAnova = pd.concat([allStats.ix[ind], dfScore0], axis=1, sort=False)

# writer = pd.ExcelWriter('test_output/ANOVA_file.xlsx', engine='xlsxwriter')
# dfAnova.to_excel(writer, sheet_name='ANOVA')

# writer.save()