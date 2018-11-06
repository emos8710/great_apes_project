import gzip
"""
NC_0014222 is phiX = not a true signal but a virus added later to the s
sample
"""

sample_file = "sample_data/Gorilla_beringei_beringei-Imfura.mpile.gz"
virus_size_file = "sample_data/virus_genome_sizes.tsv"

virus_sizes = {}
with open(virus_size_file) as f1:
    for line in f1:
        print line
        (key, value) = line.strip().split(" ")
        virus_sizes[key] = value

#with gzip.open(filename) as f2:
#    for line in f2:
#        splitted = line.strip().split("\t")
#        print splitted
