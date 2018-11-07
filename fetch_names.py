from Bio import Entrez

file = open("sample_data/virus_genome_sizes.tsv")

virus_id = []
for line in file:
    split = line.strip().split(" ")
    virus_id.append(split[0])


