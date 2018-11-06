import gzip
import os

files = []
for name in os.listdir('sample_data'):
	files.append(gzip.open('sample_data/' + name))

#path = 'C:/Users/Kristina/Documents/GitHub/great_apes_project/sample_data'

file = 'sample_data/Gorilla_beringei_graueri-Serufuli.mpile.gz'

def function_A(file):
	#Function which estimate the abundance (number of copies) of each virus for an individual.

	tup = []
	groups = {}
	output = {}

	#Open file
	with gzip.open(file) as f:
		for line in f:
			split = line.strip().split("\t")
			tup.append((split[0], split[3]))

	for x,y in tup:
		group = groups.get(x, [])
		group.append(y)
		groups[x] = group

	for x in groups:
		num_cov= len(groups[x])
		lens  = [len(s) for s in groups[x]]
		num_nuc = sum(lens)
		abundance = float(num_nuc/num_cov)
		output[x] = abundance
	print output

print function_A(file)