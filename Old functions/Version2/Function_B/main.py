import gzip
import os
from function_B_v2 import functionB
from map_percent_filter import map_percent_filter

MappedThreshold = 0.1

virus_size_file = "sample_data/virus_genome_sizes.tsv"

files = []
path = 'sample_data/'   # Directory where all files are stored
for name in os.listdir(path):   # save all files to files variable
    files.append(gzip.open(path + name))

# Read virus sizes into dict
virus_sizes = {}
with open(virus_size_file) as f1:
    for line in f1:
        (key, value) = line.strip().split(" ")
        virus_sizes[key] = int(value)
f1.close()

for i in range(len(files)):
    functionB(files[i].name)
