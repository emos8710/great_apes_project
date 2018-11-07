from Bio import Entrez
import csv

# ============================

# Define input and output file
virus_filename = 'sample_data/virus_genome_sizes.tsv'
virus_output = 'sample_data/virus_names.tsv'

virus_ids = []

# Save all 9569 virus id:s
with open(virus_filename) as f:
    for line in f:
        split = line.strip().split(" ")
        virus_ids.append(split[0])

# Fetch fasta files for each virus id from NCBI:s nuccore database
Entrez.email = 'kristina.benevides@live.se'
handle = Entrez.efetch(db='nuccore', id=virus_ids, rettype="fasta", retmode="xml")
response = Entrez.read(handle, validate = False)

# Save accession number and organism name in a tsv file
with open(virus_output,'wt') as outfile:
    for entry in response:
        tsv_writer = csv.writer(outfile, delimiter='\t')
        tsv_writer.writerow((entry.get('TSeq_accver'),entry.get('TSeq_orgname')))

# ============================