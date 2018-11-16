def functionB(sample_file, virus_ids):
	import gzip
	import numpy as np
	from collections import defaultdict
	from scipy import stats
	import matplotlib.pyplot as plt

	"""
	Function B testing either with mean comparisons or a normality test. Outputs either [mean, std] or p-value.
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
	# result = {}
	# for key in maps:
	# 	temp_mean = np.mean(maps[key])
	# 	temp_std = np.std(maps[key])
	# 	exp_mean = (min(maps[key]) + max(maps[key])) / 2
	# 	if exp_mean - temp_std <= temp_mean <= exp_mean + temp_std:
	# 		good = True
	# 	else:
	# 		good = False
	# 	if good:
	# 		result[key] = [temp_mean, temp_std]

	# Use normality test on the distribution of coverage
	result = {}
	for key in maps:
		k2, p = stats.normaltest(maps[key])
		if p > 1e-3:
			result[key] = p
		else:
			bars = np.bincount(maps[key])
			temp = [x for x in bars if x != 0]
			remove_from = int(len(temp) * 0.85)
			new_maps = [x for x in maps[key] if x < remove_from]
			k2, p = stats.normaltest(new_maps)
			if p > 1e-3:
				result[key] = p

		# Plots the HHV4 virus after trimming the last peak
		# if key == 'NC_009334.1':
		# 	plt.figure()
		# 	bars = np.bincount(new_maps)
		# 	plt.bar(range(0, len(bars)), bars)
		# 	plt.show()

	# Bar plots of the coverage for several viruses.
	nr = 25 if len(result) > 25 else len(result)
	to_plot = result.keys()

	plt.figure()
	for i in range(1, nr+1):
		plt.subplot(5, 5, i)
		bars = np.bincount(maps[to_plot[i-1]])
		plt.bar(range(0, len(bars)), bars)
	plt.show()

	return result
