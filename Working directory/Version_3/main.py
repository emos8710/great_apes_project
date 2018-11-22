import gzip
import os
from collections import defaultdict
from map_percent_filter import map_percent_filter

MappedThreshold = 0.1

# Location of virus info files
virus_size_file = "virus_data/virus_genome_sizes.tsv"
virus_name_file = "virus_data/virus_names.tsv"

# SAVE SAMPLE FILE NAMES
files = []
path = 'sample_data/'  # Directory where all files are stored
for name in os.listdir(path):  # save all files to files variable
	files.append(gzip.open(path + name))

# GO THROUGH ALL FILES
Cov_result = {}
for i in range(len(files)):
	mapped = defaultdict(int)  # used for input to map_percent_filter
	with gzip.open(files[i].name) as f1:
		for line in f1:
			(seqid, loc, nuc, nucmap) = line.strip().split("\t")
			mapped[seqid] = mapped[seqid] + 1

	print 'Analysing file ', files[i].name

	# Map percent filter with dict input
	print 'Filtering viruses with low mapping %. Threshold: ', MappedThreshold
	map_result_1 = map_percent_filter(mapped, virus_size_file, MappedThreshold)

	# Map percent filter with file input
	map_result_2 = map_percent_filter(files[i].name, virus_size_file, MappedThreshold)

	print 'With dict', map_result_1
	print 'With file', map_result_2
