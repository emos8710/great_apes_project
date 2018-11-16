def map_percent_filter(sample_file, threshold, virus_sizes):
	from collections import defaultdict
	import gzip

	"""
	Takes a dict with number of mapped nucleotides per virus nucleotide, a threshold defining the lowest mapped % that 
	you want to keep (0.0 - 1.0), and a dict with virus sizes.
	Outputs a dict with the viruses that pass the threshold and the fraction of mapped nucleotides.

	Input item example: 'NC_015050.1': [1, 1, 2, 3, 4, 5, 5, 2, 4, 5, 5, 5, 3]
	Output item example: 'NC_015050.1': 0.12
	"""

	# Count how many nucleotides are mapped for each virus in the sample
	mapped = defaultdict(int)
	with gzip.open(sample_file) as f1:
		for line in f1:
			(seqid, loc, nuc, nucmap) = line.strip().split("\t")
			mapped[seqid] = mapped[seqid] + 1
	f1.close()

	# Calculate % of mapped nucs for each virus and check with threshold
	mapped_over_th = {}
	for key in mapped:
		map_perc = float(mapped[key]) / virus_sizes[key]
		if map_perc >= threshold:
			mapped_over_th[key] = map_perc

	# Return dict of virus ids with the fraction of mapped nucleotides
	return mapped_over_th
