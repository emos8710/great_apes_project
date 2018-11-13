def functionB(sample_file):
    import gzip
    import numpy as np
    from collections import defaultdict
    import matplotlib.pyplot as plt

    """
    Output: Dict with virus codes as key, value = [%unmapped, mean, std, good]
    """

    # sample_file = "sample_data/Gorilla_beringei_beringei-Imfura.mpile.gz"
    sample_file = "sample_data/Homo_sapiens_C10.mpile.gz"
    virus_size_file = "sample_data/virus_genome_sizes.tsv"

    # Read virus sizes into dict
    virus_sizes = {}
    with open(virus_size_file) as f1:
        for line in f1:
            (key, value) = line.strip().split(" ")
            virus_sizes[key] = int(value)
    f1.close()

    # Read mapped nucleotides into a dict of lists containing the number of maps per nuc
    # Example result: 'NC_015050.1': [1, 1, 2, 3, 4, 5, 5, 2, 4, 5, 5, 5, 3]
    maps = defaultdict(list)
    with gzip.open(sample_file) as f2:
        for line in f2:
            (seqid, loc, nuc, nucmap) = line.strip().split("\t")
            maps[seqid].append(len(nucmap))
    f2.close()

    # Calculate % of unmapped nucs for each virus
    unmapped = {}
    for key in maps:
        unmapped[key] = float((virus_sizes[key] - len(maps[key]))) / virus_sizes[key]

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
        stats[key] = [unmapped[key], temp_mean, temp_std, good]