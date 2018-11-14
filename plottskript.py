
import gzip
import os
import csv
import pandas as pd
from collections import defaultdict
from functionC import functionC
from functionA import functionA
from functionB import functionB
import matplotlib.pyplot as plt

CoverageThreshold = 0.9
ErrorRateThreshold = 0.3
AbundanceThreshold = 5

files = []
path = 'sample_data/' #Directory where all files are stored
for name in os.listdir(path): #save all files to files variable
	files.append(gzip.open(path + name))

#ApeDictionary = {}

#output = 'test_output/output.tsv' #Output path
headers = ['FileName','VirusID','VirusName','Abundance','Coverage','ErrorRate']

virusFile = "virus_data/virus_names.tsv"
virus = {}
with open(virusFile) as f:
    for line in f:
        (id, name) = line.strip().split("\t")
        virus[id] = name
f.close()

Apes = {}

#This function takes an opened file with tab separated values on the form
#VirusId	Position	ReferenceNucleotide	MappedNucleotides
# and returns, for each virus in the file, the virus ID and the error rate of mapping (wrongly mapped nucleotides / all mapped nucelotides) as a dictionary
for i in range(len(files)): #go through all files and call functions A, B and C.
	viruslist = []
	viruses = {}
	with gzip.open(files[i].name) as apfil:
		for line in apfil:
			(seqid, loc, nuc, nucmap) = line.strip().split("\t") #read the lines and divide by column
			viruslist.append((seqid, nuc, nucmap, loc))
	apfil.close()

	for x,y,z,a in viruslist: #go through four elements at a time
		virus = viruses.get(x, [])
		virus.append([y, z, a])
		viruses[x]= virus

		#print len(viruses[viruses.keys()[0]])
	mismatches = {}
	errorRate ={}
	for x in viruses:
		if x=='NC_010153.1':
			counter = 0
			misnuc = 0
			Fractions=[]
			Positions=[]
			for m in range(len(viruses.get(x))): #Go through all virus postitions present in the sample (all rows in the table belonging to this virus)
				mapping = viruses.get(x)[m] #List of mapping, example ['T', 'TT', '3257'] = [refNuc, MappedNuc, position]
				maps=mapping[1] #The mapped nucleotides, 'TT'
				A=0
				T=0
				C=0
				G=0		
				for j in range(len(maps)): #go through the mapped nucleotides
					counter = counter +1
					#if maps[j] != mapping[0]: #check if mapped nucs are same as reference
					#	misnuc = misnuc +1
						#Check how coherrent the maps are

					if maps[j]=='A':
						A = A +1
					elif maps[j]=='T':
						T = T +1
					elif maps[j]=='C':
						C = C +1
					elif maps[j]=='G':
						G = G +1

				DominantNuc = max(A,T,C,G) #Find the number of the dominant nucleotide
				DisagreeingNucs = sum([A,T,C,G]) - DominantNuc #Find how many mapped nucs are not the dominant one
				Fraction = float(DisagreeingNucs)/float(sum([A,T,C,G]))  #Calcualte the fraction of nucs that are not dominant
				#print [DominantNuc, DisagreeingNucs]
				Fractions.append(Fraction)
				Positions.append(mapping[2])
			plt.bar(Positions, Fractions)
			plt.show()

					#if x=='NC_009334.1':


#plt.plot(Positions, Fractions)
#plt.show()
#print Positions
#print Fractions
		#mismatches[x] = ([counter, misnuc])
		#errorRate[x] = float(misnuc)/float(counter)
		#Apes[files[i].name] = errorRate



#NC_009334.1



