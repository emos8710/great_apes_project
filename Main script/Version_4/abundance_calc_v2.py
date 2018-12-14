def functionA(virus_maps):

	"""
	Calculates abundance.

	Parameters
	----------
	virus_maps : list
		List where each pair of items are the virus ID and the mapped nucleotides for one position.
		Example:
			('NC_001330.1', AAAA', 'NC_001330.1', 'TTTTATT')

	Returns
	-------
	dict
		Dictionary with virus ID's as keys and the abundance as values.

	"""

	# Create dictionary with virus IDs as keys and a list of the mapped nucleotides as values.
	groups = {}
	for ID, mapnuc in virus_maps:
		group = groups.get(ID, [])
		group.append(mapnuc)
		groups[ID] = group

	# Create dictionary with NCBI virus id:s as keys and abundance of virus (number of virus copies present) as values.
	# For example: {'NC_008892.1': 6.2272727272727275, 'NC_006276.1': 2.65, 'NC_016447.1': 1.5583756345177664, ...}
	output = {}
	for x in groups:
		num_cov = len(groups[x])
		lens = [len(s) for s in groups[x]]
		num_nuc = sum(lens)
		abundance = float(num_nuc)/float(num_cov)
		output[x] = abundance

	return output
