import gzip
import os
import csv
import pandas as pd
from collections import defaultdict
from map_percent_filter_v3 import map_percent_filter
from Incoherence_filter_wsize import incoherence_filter_sw
from functionA_v2 import functionA
from function_B_v3 import functionB
import numpy as np
from functionC_v2 import functionC
import matplotlib.pyplot as plt


MappedThreshold = 0.1
incoherenceThreshold = 0.05
ErrorRateThreshold = 0.03
AbundanceThreshold = 7
sliding_windowsize = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
jumpSize = 100
smooth = 1  # change to 0 if you do not want sliding window averaging of sites in the incoherence filtering.
window_score = []
zero_score = []
nonzero_score = []
delete_score = []

# Location of virus info files
virus_size_file = "virus_data/virus_genome_sizes.tsv"
virus_name_file = "virus_data/virus_names.tsv"

# SAVE SAMPLE FILE NAMES
files = []
path = 'sample_data/'  # Directory where all files are stored
for name in os.listdir(path):  # save all files to files variable
	files.append(gzip.open(path + name))

# READ VIRUS SIZES INTO DICT
virus_sizes = {}
with open(virus_size_file) as f1:
	for line in f1:
		(key, value) = line.strip().split(" ")
		virus_sizes[key] = int(value)
f1.close()

# FETCH VIRUS NAMES
virus = {}
with open(virus_name_file) as f:
	for line in f:
		(id, name) = line.strip().split("\t")
		virus[id] = name
f.close()

# GO THROUGH ALL FILES
Cov_result = {}
for i in range(len(files)):
	print 'Analysing file ', files[i].name

	# Filter out viruses with less than given mapped threshold to speed up later computations
	print 'Filtering viruses with low mapping %. Threshold: ', MappedThreshold
	MPF_result = map_percent_filter(files[i].name, virus_size_file, MappedThreshold)

	virus_ids = MPF_result.keys()  # Save the IDs of viruses that passed

	# Filter the remaining virus data based on coverage, using functionB
	print 'Trimming sites with unusually high coverage.'
	Cov_result = functionB(files[i].name, virus_ids, mean_th=1.0, peak_hicov=0.15, peak_locov=0.15)

	# Filter the data again, this time based on incoherence in the mapping
	print 'Trimming sites with incoherent mapping. Threshold: ', incoherenceThreshold
	Resultlist = []
	zeroList = []
	nonzeroList = []
	deleteList = []
	for j in sliding_windowsize:
		finalResult, zero, nonzero, delete = incoherence_filter_sw(Cov_result, incoherenceThreshold, smooth, j, jumpSize)
		Resultlist.append(finalResult)
		zeroList.append(zero)
		nonzeroList.append(nonzero)
		deleteList.append(delete)
	window_score.append(Resultlist)
	zero_score.append(zeroList)
	nonzero_score.append(nonzeroList)
	delete_score.append(deleteList)

result = np.mean(window_score, axis=0)
zero_result = np.mean(zero_score, axis=0)
zero_result = np.mean(zero_result, axis=1)
nonzero_result = np.mean(nonzero_score, axis=0)
nonzero_result = np.mean(nonzero_result, axis=1)
delete_result = np.mean(delete_score, axis=0)

f, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4)
ax1.plot(sliding_windowsize, result)
ax2.plot(sliding_windowsize, nonzero_result)
ax3.plot(sliding_windowsize, delete_result)
ax4.plot(sliding_windowsize, zero_result)

ax1.set_title('combined result')
ax2.set_title('kept nonzeros')
ax3.set_title('#removed zeros/#removed')
ax4.set_title('#removed zeros')

plt.show()

