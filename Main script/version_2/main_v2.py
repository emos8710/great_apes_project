import gzip
import os
import csv
import pandas as pd
from collections import defaultdict
from map_percent_filter_v2 import map_percent_filter
from Incoherence_filter import Incoherence_filter
from functionA_v2 import functionA
from function_B_v3 import functionB
from functionC_v2 import functionC


MappedThreshold = 0.1
incoherenceThreshold = 0.07
incoherenceThresholdSmooth = 0.03
CoverageThreshold = 0.9
ErrorRateThreshold = 0.3
AbundanceThreshold = 10
binsize = 3
jumpSize = 100

output = 'test_output/output.tsv'		# Output path
headers = ['FileName', 'VirusID', 'VirusName', 'Abundance', 'Coverage', 'ErrorRate']  # headers in the file

virus_size_file = "virus_data/virus_genome_sizes.tsv"  # where to find file with virus genome sizes

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
virusFile = "virus_data/virus_names.tsv"
virus = {}
with open(virusFile) as f:
	for line in f:
		(id, name) = line.strip().split("\t")
		virus[id] = name
f.close()

# GO THROUGH ALL FILES
Cov_result = {}
for i in range(len(files)):
	mapped = defaultdict(int)  # used for input to map:percent_filter
	with gzip.open(files[i].name) as f1:
		for line in f1:
			(seqid, loc, nuc, nucmap) = line.strip().split("\t")
			mapped[seqid] = mapped[seqid] + 1
	f1.close()

	print 'Analysing file ', files[i].name

	# Filter out viruses with less than given mapped threshold to speed up later computations
	print 'Filtering viruses with low mapping %. Threshold: ', MappedThreshold
	Cov_result = map_percent_filter(mapped, MappedThreshold, virus_sizes)

	virus_ids = Cov_result.keys()  # save the IDs of viruses that passed

	# Filter the data based on coverage, using functionB
	print 'Trimming sites with unusually high coverage.'
	Cov_result = functionB(files[i].name, virus_ids, mean_th=1.0, peak_hicov=0.15, peak_locov=0.15)

	# Filter the data again, this time based on incoherence in the mapping
	print 'Trimming sites with incoherent mapping. Threshold: ', incoherenceThreshold, 'Smoothed threshold: ',\
		incoherenceThresholdSmooth

	inc_filtered = Incoherence_filter(Cov_result, incoherenceThreshold, incoherenceThresholdSmooth, binsize, jumpSize)
	# Prepare output from filtering to be input to the functions A, B and C
	inputA = []
	inputB = {}
	inputC = []
	for key in inc_filtered:
		inputB[key] = len(inc_filtered[key]['mapped_nucs'])  # change to 'mapped_nucs_smooth' for smoothed version
		for j in range(len(inc_filtered[key]['mapped_nucs'])):  # change to 'mapped_nucs_smooth' for smoothed version
			inputA.append((key, inc_filtered[key]['mapped_nucs'][j])) # change to 'mapped_nucs_smooth' for smoothed version
			inputC.append((key, inc_filtered[key]['ref_nuc'][j], inc_filtered[key]['mapped_nucs'][j]))  # change to 'ref_nuc_smooth, 'mapped_nucs_smooth' for smoothed version

print 'Filtering and printing to file. Abundance th: ', AbundanceThreshold, 'Coverage th: ', CoverageThreshold, \
	'Error rate th: ', ErrorRateThreshold

with open(output, 'wt') as outfile:
	tsv_writer = csv.writer(outfile, delimiter='\t')
	tsv_writer.writerow(headers)
	for i in range(len(files)):  # go through all files and call functions A, B and C.
		Abundance = functionA(inputA)
		Coverage = map_percent_filter(inputB, 0, virus_sizes)
		ErrorRate = functionC(inputC)
		for key in Abundance:  # Go through all virus IDs, should be same from all functions from the same file
			if Abundance.get(key) > AbundanceThreshold and Coverage.get(key) < CoverageThreshold and \
					ErrorRate.get(key) < ErrorRateThreshold:
				for id in virus:
					if id == key:
						tsv_writer.writerow((files[i].name, key, virus[id], float(round(Abundance.get(key), 2)),
											 float(round(1 - Coverage.get(key), 3)),
											 float(round(ErrorRate.get(key), 3))))
outfile.close()

niceOutput = pd.read_csv('test_output/output.tsv', sep='\t')
niceOutput.to_csv("test_output/output_excel_file.xls", sep='\t', index=False)


# final_result[files[i].name] = inc_filtered
