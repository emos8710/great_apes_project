def map_percent_filter(virus_data, virus_size_file, threshold=0.0, removed_sites=None):

	"""
	Calculates how much of the reference virus has been mapped to.

	Parameters
	----------
	virus_data : str or dict
		String with the virus data file location or dictionary of the number of mapped nucleotides per virus.
		Dict example:
			{'NC_015050.1': 500, 'NC_009127.1': 600}
	virus_size_file : str
		Location of tsv file containing virus sizes.
	threshold : float, optional
		Minimum mapped threshold. Must be in range 0.0 - 1.0.
		Value of 0.0 returns all viruses (default) and 1.0 returns no viruses.
	removed_sites : dict
		Number of sites removed in previous filtering steps for each virus. Virus ID as key, nr of sites as value.

	Returns
	-------
	dict
		Dictionary with the mapped fraction for each virus above the threshold.
		Example:
			{'NC_015050.1': 0.12, 'NC_009127.1': 0.30}

	Raises
	------
	ValueError
		If virus_data is not str or dict.
	ValueError
		If threshold is not is span [0.0, 1.0]
	"""

	import gzip
	from collections import defaultdict

	# Check that threshold is valid
	if not 0.0 <= threshold <= 1.0:
		raise ValueError('threshold must be a value in the span [0.0, 1.0].')

	# Read virus sizes from the file
	virus_sizes = {}
	with open(virus_size_file) as f1:
		for line in f1:
			(seqid, size) = line.strip().split(' ')
			virus_sizes[seqid] = int(size)

	# Remove removed sites from virus lengths
	if isinstance(removed_sites, dict):
		for key in virus_sizes:
			virus_sizes[key] = virus_sizes[key] - removed_sites[key]

	mapped_over_th = {}

	# Calculate mapped fraction if data is in dict
	if isinstance(virus_data, dict):
		for key in virus_data:
			map_perc = float(virus_data[key]) / virus_sizes[key]
			if map_perc >= threshold:
				mapped_over_th[key] = map_perc

	# Read file and then calculate mapped fraction if data is in file
	elif isinstance(virus_data, str):
		mapped = defaultdict(int)  # used for input to map_percent_filter
		with gzip.open(virus_data) as f2:
			for line in f2:
				(seqid, loc, nuc, nucmap) = line.strip().split('\t')
				mapped[seqid] = mapped[seqid] + 1

		for key in mapped:
			map_perc = float(mapped[key]) / virus_sizes[key]
			if map_perc >= threshold:
				mapped_over_th[key] = map_perc

	# Raise error if virus_data is neither dict nor string
	else:
		raise ValueError('virus_data must be dict or str.')

	# Return dict of virus ids with the fraction of mapped nucleotides
	return mapped_over_th
