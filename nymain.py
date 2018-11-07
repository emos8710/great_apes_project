import gzip
import os
from collections import defaultdict
from functionC import functionC
from functionA import functionA
from functionB import functionB

files = []
path = 'sample_data/' #Directory where all files are stored
for name in os.listdir(path): #save all files to files variable
	files.append(gzip.open(path + name))

for i in range(len(files)): #go through all files and call functions A, B and C.
	Abundance = functionA(files[i].name)
	Coverage = functionB(files[i].name)
	ErrorRate = functionC(files[i].name)
	files[i]name.close() #Close the file
	


coverageThreshold = 50
binomialThreshold = 
ErrorRateThreshold = 0.5
AbundanceThreshold = 5






