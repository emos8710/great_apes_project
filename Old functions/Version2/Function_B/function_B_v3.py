def functionB(sample_file, virus_ids):
	import gzip
	import numpy as np
	from collections import defaultdict
	import matplotlib.pyplot as plt

	"""
	Removes maps with unusually high coverage.
	"""

	# Read mapped nucleotides into a dict of lists containing the number of maps per nuc
	# Example result: 'NC_015050.1': [1, 1, 2, 3, 4, 5, 5, 2, 4, 5, 5, 5, 3]
	maps = defaultdict(list)
	with gzip.open(sample_file) as f1:
		for line in f1:
			(seqid, loc, nuc, nucmap) = line.strip().split("\t")
			if seqid in virus_ids:
				maps[seqid].append(len(nucmap))
	f1.close()

	trimmed_maps = {}
	trimmed = {}
	for key in maps:
		bars = np.bincount(maps[key])
		temp = [x for x in bars if x != 0]
		remove_from = int(len(temp) * 0.85)
		remove_to = len(temp) - 1
		mean_start = np.mean(temp[0:remove_from-1])
		if remove_from != remove_to:
			mean_end = np.mean(temp[remove_from:remove_to])
		else:
			mean_end = 0
		if mean_end > mean_start:
			new_maps = [x for x in maps[key] if x < remove_from]
			trimmed_maps[key] = new_maps
			trimmed[key] = maps[key]
		else:
			trimmed_maps[key] = maps[key]

	# Bar plots of the coverage for several viruses before trimming.
	# nr = 25 if len(maps.keys()) > 25 else len(maps.keys())
	# to_plot = maps.keys()
	#
	# plt.figure()
	# for i in range(1, nr+1):
	# 	plt.subplot(5, 5, i)
	# 	bars = np.bincount(maps[to_plot[i-1]])
	# 	plt.bar(range(0, len(bars)), bars)
	# plt.show(block=False)
	#
	# # Bar plots of the coverage for several viruses after trimming.
	# plt.figure()
	# for i in range(1, nr + 1):
	# 	plt.subplot(5, 5, i)
	# 	bars = np.bincount(trimmed_maps[to_plot[i - 1]])
	# 	plt.bar(range(0, len(bars)), bars)
	# plt.show()

	# Bar plots of the trimmed viruses pre timming
	nr = 9 if len(trimmed.keys()) > 9 else len(trimmed.keys())
	to_plot = trimmed.keys()
	plt.figure()
	for i in range(1, nr + 1):
		plt.subplot(3, 3, i)
		bars = np.bincount(trimmed[to_plot[i - 1]])
		plt.bar(range(0, len(bars)), bars)
	plt.show()

	return trimmed_maps
