def map_percent_filter(mapped, threshold, virus_sizes):
	from collections import defaultdict

	"""
	Takes a dict with number of mapped nucleotides per virus nucleotide, a threshold defining the lowest mapped % that 
	you want to keep (0.0 - 1.0), and a dict with virus sizes.
	Outputs a dict with the viruses that pass the threshold and the fraction of mapped nucleotides.

	Input item example: 'NC_015050.1': [1, 1, 2, 3, 4, 5, 5, 2, 4, 5, 5, 5, 3]
	Output item example: 'NC_015050.1': 0.12
	"""

	# Calculate % of mapped nucs for each virus and check with threshold
	mapped_over_th = {}
	for key in mapped:
		map_perc = float(mapped[key]) / virus_sizes[key]
		if map_perc >= threshold:
			mapped_over_th[key] = map_perc

	# Return dict of virus ids with the fraction of mapped nucleotides
	return mapped_over_th
