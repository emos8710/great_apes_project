def functionC(sample_file):

#This function takes an opened file with tab separated values on the form
#VirusId	Position	ReferenceNucleotide	MappedNucleotides
# and returns, for each virus in the file, the virus ID and the error rate of mapping (wrongly mapped nucleotides / all mapped nucelotides) as a dictionary
	import gzip
	viruslist = []
	viruses = {}
	with gzip.open(sample_file) as apfil:
		for line in apfil:
			(seqid, loc, nuc, nucmap) = line.strip().split("\t") #read the lines and divide by column
			viruslist.append((seqid, nuc, nucmap))
	apfil.close()

	for x,y,z in viruslist: #go through three elements at a time
		virus = viruses.get(x, [])
		virus.append([y, z])
		viruses[x]= virus

	#print len(viruses[viruses.keys()[0]])
	mismatches = {}
	errorRate ={}
	for x in viruses:
		counter = 0
		misnuc = 0
		for i in range(len(viruses.get(x))): #Go through all virus postitions present in the sample (all rows in the table belonging to this virus)
			mapping = viruses.get(x)[i] #List of mapping, example ['T', 'TT']
			maps=mapping[1] #The mapped nucleotides, 'TT'
			for j in range(len(maps)): #go through the mapped nucleotides
				counter = counter +1
				if maps[j] != mapping[0]: #check if mapped nucs are same as reference
					misnuc = misnuc +1
		mismatches[x] = ([counter, misnuc])
		errorRate[x] = [float(misnuc)/float(counter)]


	return errorRate


