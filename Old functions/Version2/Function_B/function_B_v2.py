def functionB(sample_file, virus_ids):
	import gzip
	import numpy as np
	from collections import defaultdict

	"""
	Output: Dict with virus codes as key, value = [%unmapped, mean, std, good]
	"""

	# Read mapped nucleotides into a dict of lists containing the number of maps per nuc
	# Example result: 'NC_015050.1': [1, 1, 2, 3, 4, 5, 5, 2, 4, 5, 5, 5, 3]
	maps = defaultdict(list)
	with gzip.open(sample_file) as f2:
		for line in f2:
			(seqid, loc, nuc, nucmap) = line.strip().split("\t")
			maps[seqid].append(len(nucmap))
	f2.close()

	# Calculate mean, standard deviation and expected mean of data.
	# Data with mean within one std of the expected mean are good.
	stats = {}
	for key in maps:
		temp_mean = np.mean(maps[key])
		temp_std = np.std(maps[key])
		exp_mean = (min(maps[key]) + max(maps[key])) / 2
		if (temp_mean <= exp_mean + temp_std) and (temp_mean >= exp_mean - temp_std):
			good = 1
		else:
			good = 0
		stats[key] = [temp_mean, temp_std, good]
