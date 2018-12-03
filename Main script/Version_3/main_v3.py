import gzip
import os
import csv
import pandas as pd
from collections import defaultdict
from map_percent_filter_v3 import map_percent_filter
from Incoherence_filter_SW import incoherence_filter_sw
from functionA_v2 import functionA
from function_B_v3 import functionB
from functionC_v2 import functionC


MappedThreshold = 0.1
incoherenceThreshold = 0.03
ErrorRateThreshold = 0.03
AbundanceThreshold = 7
sliding_windowsize = 50
jumpSize = 100
smooth = 1  # change to 0 if you do not want sliding window averaging of sites in the incoherence filtering.

# Create output file. OVERWRITES PREVIOUS OUTPUT FILE!
outfile_path = 'test_output/output.tsv'		# Output path
headers = ['FileName', 'VirusID', 'VirusName', 'Abundance', 'Coverage', 'ErrorRate']  # headers in the file
with open(outfile_path, 'w') as outfile:
	tsv_writer = csv.writer(outfile, delimiter='\t')
	tsv_writer.writerow(headers)

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

	inc_filtered = incoherence_filter_sw(Cov_result, incoherenceThreshold, smooth, sliding_windowsize, jumpSize)


	# Prepare output from filtering to be input to the functions A, B and C
	inputA = []
	inputB = {}
	inputC = []
	for key in inc_filtered:
		inputB[key] = len(inc_filtered[key]['mapped_nucs'])
		for j in range(len(inc_filtered[key]['mapped_nucs'])):
			inputA.append((key, inc_filtered[key]['mapped_nucs'][j]))
			inputC.append((key, inc_filtered[key]['ref_nuc'][j], inc_filtered[key]['mapped_nucs'][j]))

	print 'Filtering and printing to file. Abundance th: ', AbundanceThreshold, 'Mapped th: ', MappedThreshold, \
		'Error rate th: ', ErrorRateThreshold

	Abundance = functionA(inputA)
	Mapped = map_percent_filter(inputB, virus_size_file)
	ErrorRate = functionC(inputC)

	with open(outfile_path, 'a') as outfile:
		tsv_writer = csv.writer(outfile, delimiter='\t')
		for key in Abundance:  # Go through all virus IDs, should be same from all functions from the same file
			if Abundance.get(key) > AbundanceThreshold and Mapped.get(key) > MappedThreshold and \
					ErrorRate.get(key) < ErrorRateThreshold:
				tsv_writer.writerow([files[i].name, key, virus[key],
									float(round(Abundance.get(key), 2)),
									float(round(Mapped.get(key), 3)),
									float(round(ErrorRate.get(key), 3))])

niceOutput = pd.read_csv('test_output/output.tsv', sep='\t')
niceOutput.to_csv("test_output/output_excel_file.xls", sep='\t', index=False)

