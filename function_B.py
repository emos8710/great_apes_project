import gzip
from collections import defaultdict

"""
NC_0014222 is phiX = not a true signal but a virus added later to the s
sample
"""

sample_file = "sample_data/Gorilla_beringei_beringei-Imfura.mpile.gz"
virus_size_file = "sample_data/virus_genome_sizes.tsv"

# Read virus sizes into dict
virus_sizes = {}
with open(virus_size_file) as f1:
    for line in f1:
        (key, value) = line.strip().split(" ")
        virus_sizes[key] = int(value)
f1.close()

# Read mapped nucleotides into a dict of lists containing the number of maps per nuc
# Example result: 'NC_015050.1': [1, 1, 2, 3, 4, 5, 5, 2, 4, 5, 5, 5, 3]
maps = defaultdict(list)
with gzip.open(sample_file) as f2:
    for line in f2:
        (seqid, loc, nuc, nucmap) = line.strip().split("\t")
        maps[seqid].append(len(nucmap))

# Calculate nr of unmapped nucs for each virus
unmapped = {}
for key in maps:
    unmapped[key] = virus_sizes[key] - len(maps[key])

