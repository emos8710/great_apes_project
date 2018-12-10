def functionB(sample_file, virus_ids, mean_th=1.0, peak_locov=0.15, peak_hicov=0.15):

	"""
	Removes mapped sites with unusually high coverage.

	Example plot with peak that should be removed:
	# sites
		8|
		7|                       |
		6| |                     ||
		5| ||                   |||
		4| ||                   |||
		3| |||                 |||||
		2| |||||||||||||||||||||||||
		1|___________________________
		  0 1 2 3 4 5 6 7 8 9 10 11		coverage

		  <--->               <--->
		peak_locov            peak_hicov

	Parameters
	----------
	sample_file : str
		Location of the virus data file.
	virus_ids : list
		List containing virus ids (accession nrs) for the viruses that should be analysed. Viruses not in this list will
		be excluded from the output and the analysis.
	mean_th : float, optional
		Sets multiplier for the mean for determining threshold. If mean(peak) is greater than mean(rest) * mean_th
		the peak is removed from the data.
		mean_th > 1 means the peak must be taller before it is removed.
		mean_th < 1 means shorter peaks are removed.
	peak_locov : float, optional
		How many % of the coverage interval a low coverage peak is expected to cover.
	peak_hicov : float, optional
		How many % of the coverage interval a high coverage peak is expected to cover.

	Returns
	-------
	dict
		Dictionary with virus ID's as keys. The values are one dict per virus containing the keys trim_map_loc,
		trim_ref_nuc, trim_map_nuc and nr_removed_sites. The values of these dicts are lists with the mapping locations,
		the reference nucleotide, the mapped nucleotides and the number of sites that have been removed respectively.
	"""

	import gzip
	import numpy as np
	from collections import defaultdict
	# import matplotlib.pyplot as plt

	# Read mapped nucleotides into a dict of lists containing the number of maps per nuc
	# As well as saving the locations of each map in the virus genome
	# Example result:
	# Maps: 'NC_015050.1': [1, 1, 2, 3, 4, 5, 5, 2, 4, 5, 5, 5, 3]
	# Locations: 'NC_015050.1': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
	maps = defaultdict(list)
	map_loc = defaultdict(list)
	ref_nuc = defaultdict(list)
	map_nuc = defaultdict(list)
	with gzip.open(sample_file) as f1:
		for line in f1:
			(seqid, loc, nuc, nucmap) = line.strip().split("\t")
			if seqid in virus_ids:
				maps[seqid].append(len(nucmap))
				map_loc[seqid].append(int(loc))
				ref_nuc[seqid].append(nuc)
				map_nuc[seqid].append(nucmap)
	f1.close()

	virus_result = {}
	for key in maps:
		virus_result[key] = {}
		bars = np.bincount(maps[key])		# bars = nr of sites that have the same coverage
		temp = [x for x in enumerate(bars) if x[1] != 0]		# temp = [coverage, nr of sites]

		# Find the index of the right edge of each interval
		lim15 = int(len(temp) * peak_locov)		# lim = index marking 15%, 85% and 100% of temp
		lim85 = int(len(temp) * (1-peak_hicov))
		lim100 = len(temp) - 1

		# Calculate mean nr of sites in each coverage interval
		mean_beg = np.mean(temp[0:lim15], axis=0)[1] if lim15 != 0 else 0		# mean = [mean nr of sites in interval]
		mean_mid = np.mean(temp[lim15+1:lim85], axis=0)[1] if lim15+1 < lim85 else 0
		mean_end = np.mean(temp[lim85+1:lim100], axis=0)[1] if lim85+1 < lim100 else 0

		# Remove sites with too high coverage if there is a peak
		if mean_end > mean_mid * mean_th and mean_beg > mean_mid * 1.2:		# The 1.2 prevents scenarios where a low-cov
																			# peak is identified by mistake
			inx = [i for i, x in enumerate(maps[key]) if x < temp[lim85][0]]		# x < temp[lim85][0] removes the
																					# 15 % most covered sites
			removed = [nr for cov, nr in temp if cov >= temp[lim85][0]]
		else:
			inx = [i for i, x in enumerate(maps[key])]
			removed = [0]

		virus_result[key]['trim_map_loc'] = [map_loc[key][i] for i in inx]
		virus_result[key]['trim_ref_nuc'] = [ref_nuc[key][i] for i in inx]
		virus_result[key]['trim_map_nuc'] = [map_nuc[key][i] for i in inx]
		virus_result[key]['nr_removed_sites'] = sum(removed)

	# Bar plots of the coverage for several viruses before trimming.
	# nr = 15 if len(trimmed.keys()) > 15 else len(trimmed.keys())
	# to_plot = trimmed.keys()
	#
	# plt.figure()
	# for i in range(1, nr+1):
	# 	plt.subplot(5, 3, i)
	# 	bars = np.bincount(maps[to_plot[i-1]])
	# 	plt.bar(range(0, len(bars)), bars)
	# plt.show(block=False)

	# Bar plots of the coverage for several viruses after trimming.
	# plt.figure()
	# for i in range(1, nr + 1):
	# 	plt.subplot(5, 3, i)
	# 	bars = np.bincount(trimmed_maps[to_plot[i - 1]])
	# 	plt.bar(range(0, len(bars)), bars)
	# plt.show()

	# Bar plots of the trimmed viruses pre timming
	# nr = 9 if len(trimmed.keys()) > 9 else len(trimmed.keys())
	# to_plot = trimmed.keys()
	# plt.figure()
	# for i in range(1, nr + 1):
	# 	plt.subplot(3, 3, i)
	# 	bars = np.bincount(trimmed[to_plot[i - 1]])
	# 	plt.bar(range(0, len(bars)), bars)
	# plt.show()

	return virus_result
