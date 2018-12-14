def functionC(viruslist):

	"""
	Calculates the mapping error rate. The error rate is calculated as:
	(wrongly mapped nucleotides / all mapped nucelotides)

	Parameters
	----------
	viruslist : list
		A list where each three consecutive items are the virus ID, reference nucleotide and mapped nucleotides.
		Example:
			('NC_001330.1', 'A', 'AAAA', 'NC_001330.1', 'T', 'TTTTATT')

	Returns
	-------
	dict
		Dictionary with virus ID's as keys. The values are the error rate.
	"""

	# Create dict with virus IDs as keys and list of [refnuc, mapnuc] as values
	viruses = {}
	for ID, refnuc, mapnuc in viruslist:		# go through three elements at a time
		virus = viruses.get(ID, [])
		virus.append([refnuc, mapnuc])
		viruses[ID] = virus

	# Calculate error rate for each virus
	mismatches = {}
	errorRate = {}
	for x in viruses:
		counter = 0
		misnuc = 0
		for i in range(len(viruses.get(x))):		# Go through all virus postitions present in the sample
													# (all rows in the table belonging to this virus)
			mapping = viruses.get(x)[i]		# List of mapping, example ['T', 'TT']
			maps = mapping[1]		# The mapped nucleotides, 'TT'
			for j in range(len(maps)):		# go through the mapped nucleotides
				counter = counter + 1
				if maps[j] != mapping[0]:		# check if mapped nucs are same as reference
					misnuc = misnuc + 1
		mismatches[x] = ([counter, misnuc])
		errorRate[x] = float(misnuc)/float(counter)

	return errorRate
