import gzip
"""
NC_0014222 is phiX = not a true signal but a virus added later to the s
sample
"""

filename = "sample_data/Gorilla_beringei_beringei-Imfura.mpile.gz"

with gzip.open(filename) as f1:
    for line in f1:
        splitted = line.strip().split("\t")
        print splitted
