import gzip
import os
import csv
from map_percent_filter_v4 import map_percent_filter
from function_B_v4 import functionB

MappedThreshold = 0.1

# Location of virus info files
virus_size_file = "virus_data/virus_genome_sizes.tsv"
virus_name_file = "virus_data/virus_names.tsv"

# SAVE SAMPLE FILE NAMES
files = []
path = 'test_sample/'  # Directory where all files are stored
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

	removed_sites = {}
	for key in Cov_result:
		removed_sites[key] = Cov_result[key]['nr_removed_sites']

	input_MPF = {}
	for key in Cov_result:
		input_MPF[key] = len(Cov_result[key]['trim_map_nuc'])

	MPF_result_2 = map_percent_filter(input_MPF, virus_size_file, removed_sites=removed_sites)
	MPF_result_3 = map_percent_filter(input_MPF, virus_size_file)

	compare = {}
	for key in MPF_result_2:
		if MPF_result_2[key] != MPF_result_3[key]:
			compare[key] = [MPF_result_2[key], MPF_result_3[key]]

	print compare
