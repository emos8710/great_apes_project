def functionB(sample_file, virus_ids, mean_th=1.0, peak_locov=0.15, peak_hicov=0.15):
	import gzip
	import numpy as np
	from collections import defaultdict
	import matplotlib.pyplot as plt

	"""
	Removes maps with unusually high coverage.
	mean_th sets how much above the mean the last peak must be to be removed.
	peak_locov sets how much of the coverage interval the low coverage peak is expected to cover.
	peak_hicov sets how much of the coverage interval the high coverage peak is expected to cover.
	"""

	# Read mapped nucleotides into a dict of lists containing the number of maps per nuc
	# As well as saving the locations of each map in the virus genome
	# Example result:
	# Maps: 'NC_015050.1': [1, 1, 2, 3, 4, 5, 5, 2, 4, 5, 5, 5, 3]
	# Locations: 'NC_015050.1': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
	maps = defaultdict(list)
	map_locs = defaultdict(list)
	with gzip.open(sample_file) as f1:
		for line in f1:
			(seqid, loc, nuc, nucmap) = line.strip().split("\t")
			if seqid in virus_ids:
				maps[seqid].append(len(nucmap))
				map_locs[seqid].append(loc)

	f1.close()

	trimmed_maps = {}
	trimmed = {}
	for key in maps:
		bars = np.bincount(maps[key])		# bars = [nr of sites]
		temp = [x for x in enumerate(bars, 0) if x[1] != 0]		# temp = [coverage, nr of sites]
		lim15 = int(len(temp) * peak_locov)		# lim = index marking 15%, 85% and 100% of temp
		lim85 = int(len(temp) * (1-peak_hicov))
		lim100 = len(temp) - 1
		mean_beg = np.mean(temp[0:lim15], axis=0)[1] if lim15 != 0 else 0		# mean = [mean nr of sites in interval]
		mean_mid = np.mean(temp[lim15+1:lim85], axis=0)[1] if lim15+1 < lim85 else 0
		mean_end = np.mean(temp[lim85+1:lim100], axis=0)[1] if lim85+1 < lim100 else 0
		if mean_end > mean_mid * mean_th and mean_beg > mean_mid:
			new_maps = [x for x in maps[key] if x < temp[lim85][0]]		# x < temp[lim85][0] removes the 15 % most covered sites
			trimmed_maps[key] = new_maps
			trimmed[key] = maps[key]
		else:
			trimmed_maps[key] = maps[key]

	print len(trimmed)

	# Bar plots of the coverage for several viruses before trimming.
	nr = 15 if len(trimmed.keys()) > 15 else len(trimmed.keys())
	to_plot = trimmed.keys()

	plt.figure()
	for i in range(1, nr+1):
		plt.subplot(5, 3, i)
		bars = np.bincount(maps[to_plot[i-1]])
		plt.bar(range(0, len(bars)), bars)
	plt.show(block=False)

	# Bar plots of the coverage for several viruses after trimming.
	plt.figure()
	for i in range(1, nr + 1):
		plt.subplot(5, 3, i)
		bars = np.bincount(trimmed_maps[to_plot[i - 1]])
		plt.bar(range(0, len(bars)), bars)
	plt.show()

	# Bar plots of the trimmed viruses pre timming
	# nr = 9 if len(trimmed.keys()) > 9 else len(trimmed.keys())
	# to_plot = trimmed.keys()
	# plt.figure()
	# for i in range(1, nr + 1):
	# 	plt.subplot(3, 3, i)
	# 	bars = np.bincount(trimmed[to_plot[i - 1]])
	# 	plt.bar(range(0, len(bars)), bars)
	# plt.show()

	return trimmed_maps
