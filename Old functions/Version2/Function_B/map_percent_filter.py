def map_percent_filter(maps, threshold, virus_sizes):
    """
    Takes a dict with number of mapped nucleotides per virus nucleotide, a threshold (0 - 1) and  dict with virus sizes.
    Outputs a dict with the viruses that pass the threshold and their mapped %.

    Input item example: 'NC_015050.1': [1, 1, 2, 3, 4, 5, 5, 2, 4, 5, 5, 5, 3]
    Output item example: 'NC_015050.1': 0.99
    """

    # Calculate % of mapped nucs for each virus and check with threshold
    mapped_over_th = {}
    for key in maps:
        mapped = float(len(maps[key])) / virus_sizes[key]
        if mapped >= threshold:
            mapped_over_th[key] = mapped

    return mapped_over_th
