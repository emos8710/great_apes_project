def functionA(file):
	import gzip
	# Function which estimate the abundance (number of copies) of each virus for an individual and stores this information
	# in a dictionary.

	tup = []
	groups = {}
	output = {}

	# Assumes file format with rows containing NCBI virus id, position on the virus genome, virus genome nucelotide,
	# nucleotides found in the great ape sample, for example ['NC_001993.1', '229479', 'T', 'TTTT'].
	with gzip.open(file) as f:
		for line in f:
			split = line.strip().split("\t")
			tup.append((split[0], split[3]))
	f.close()		
    # Dictionary with NCBI virus id:s as keys and abundance of virus (number of virus copies present) as values.
	# For example: {'NC_008892.1': 6.2272727272727275, 'NC_006276.1': 2.65, 'NC_016447.1': 1.5583756345177664, ...}
	for x,y in tup:
		group = groups.get(x, [])
		group.append(y)
		groups[x] = group

	for x in groups:
		num_cov= len(groups[x])
		lens  = [len(s) for s in groups[x]]
		num_nuc = sum(lens)
		abundance = float(num_nuc)/float(num_cov)
		output[x] = abundance
	return output

