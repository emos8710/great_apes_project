def functionC(viruslist):

# This function takes a list where each three consecutive positions contain
# (VirusId,	ReferenceNucleotide, 	MappedNucleotides)
# and returns, for each virus in the list, the virus ID and the error rate of mapping (wrongly mapped nucleotides / all mapped nucelotides) as a dictionary
	viruses = {}

	for x,y,z in viruslist: #go through three elements at a time
		virus = viruses.get(x, [])
		virus.append([y, z])
		viruses[x] = virus

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
		errorRate[x] = float(misnuc)/float(counter)


	return errorRate


