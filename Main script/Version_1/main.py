import gzip
import os
import csv
import pandas as pd
from collections import defaultdict
from functionC import functionC
from functionA import functionA
from functionB import functionB

CoverageThreshold = 0.9
ErrorRateThreshold = 0.3
AbundanceThreshold = 5

files = []
path = '../../sample_data/'  # Directory where all files are stored
for name in os.listdir(path):  # Save all files to files variable
	files.append(gzip.open(path + name))

# ApeDictionary = {}

output = '../../test_output/output.tsv'  # Output path
headers = ['FileName', 'VirusID', 'VirusName', 'Abundance', 'Coverage', 'ErrorRate']

virusFile = "../../virus_data/virus_names.tsv"
virus = {}
with open(virusFile) as f:
    for line in f:
        (id, name) = line.strip().split("\t")
        virus[id] = name
f.close()

with open(output, 'wt') as outfile:
	tsv_writer = csv.writer(outfile, delimiter='\t')
	tsv_writer.writerow(headers)
	for i in range(len(files)): #go through all files and call functions A, B and C.
		Abundance = functionA(files[i].name)
		Coverage = functionB(files[i].name)
		ErrorRate = functionC(files[i].name)
		#AmazingViruses = {} #These are the ones we keep after filtering!
		for key in Abundance: #Go through all virus IDs, should be same from all functions from the same file
			#keptViruses = {} #To store the info about the Amazing viruses
			if  (Abundance.get(key)>AbundanceThreshold and Coverage.get(key)[3]==1 and Coverage.get(key)[0]<CoverageThreshold) and ErrorRate.get(key)[0]<ErrorRateThreshold:
				#keptViruses = {'Abundance': Abundance.get(key), 'Coverage': 1 - Coverage.get(key)[0], 'ErrorRate': ErrorRate.get(key)}
				#AmazingViruses[key] = keptViruses
				#ApeDictionary[files[i].name] = AmazingViruses
				for id in virus:
					if id == key:
						tsv_writer.writerow((files[i].name,key,virus[id],float(round(Abundance.get(key),2)), float(round(1-Coverage.get(key)[0],3)), float(round(ErrorRate.get(key)[0],3))))
outfile.close()

niceOutput = pd.read_csv('../../test_output/output.tsv', sep='\t')
niceOutput.to_csv("../../test_output/output_excel_file.xls", sep='\t', index=False)


#print ApeDictionary

#Ape dictionary has file names as keys and values are dictionaries with virus IDs as keys. This inner dictionary has a dictionary containing the attributes of the virus. Example with two ape files whre the first one contains one virus and the second one conatins two viruses:
#{'sample_data/testapa.mpile.gz': {'NC_014069.1': {'Abundance': 3.8292282430213467, 'ErrorRate': [0.08404802744425385], 'Coverage': 0.3300813008130081}}, 'sample_data/Gorilla_beringei_beringei-Umurimo.mpile.gz': {'NC_008913.1': {'Abundance': 5.430722891566265, 'ErrorRate': [0.04991680532445923], 'Coverage': 0.10516312955337348}, 'NC_009889.1': {'Abundance': 12.477586206896552, 'ErrorRate': [0.1494058311455023], 'Coverage': 0.2759933380918391}}}


