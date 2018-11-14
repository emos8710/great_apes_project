def functionB(sample_file, virus_ids):
	import gzip
	import numpy as np
	from collections import defaultdict
	import matplotlib.pyplot as plt

	"""
	Output: Dict with virus codes as key, value = [%unmapped, mean, std, good]
	"""

	# Read mapped nucleotides into a dict of lists containing the number of maps per nuc
	# Example result: 'NC_015050.1': [1, 1, 2, 3, 4, 5, 5, 2, 4, 5, 5, 5, 3]
	maps = defaultdict(list)
	with gzip.open(sample_file) as f2:
		for line in f2:
			(seqid, loc, nuc, nucmap) = line.strip().split("\t")
			if seqid in virus_ids:
				maps[seqid].append(len(nucmap))
	f2.close()

	# Calculate mean, standard deviation and expected mean of data.
	# Data with mean within one std of the expected mean are good.
	stats = {}
	for key in maps:
		temp_mean = np.mean(maps[key])
		temp_std = np.std(maps[key])
		exp_mean = (min(maps[key]) + max(maps[key])) / 2
		if exp_mean - temp_std <= temp_mean <= exp_mean + temp_std:
			good = True
		else:
			good = False
		if good:
			stats[key] = [temp_mean, temp_std]

	# Bar plots of the coverage for several viruses.
	nr = 25 if len(stats) > 25 else len(stats)
	to_plot = stats.keys()

	plt.figure()
	for i in range(1, nr+1):
		plt.subplot(5, 5, i)
		bars = np.bincount(maps[to_plot[i-1]])
		plt.bar(range(0, len(bars)), bars)
	plt.show()

	return stats
