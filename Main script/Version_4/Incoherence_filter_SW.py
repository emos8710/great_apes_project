def incoherence_filter_sw(Cov_result, incoherenceThreshold, smooth, binsize=50, jumpSize=100):
    """
    Removes positions or regions with high incoherence in the mapping. If smoothing is specified, a sliding window is used.
    the window will not slide across gaps in the mapping.

    Ezample of sliding window. The dotted vertical line represents a large gap in the mapping. The numbers are the
    incoherence score of those positions, and the brackets show the sliding sindows that would be applied to this virus.
                 |
    0   0   0.3  |  0.3   0   0   0.5   0       Score:
    |_________|  |                                 0.1
                    |________|                     0.1
                          |________|               0.17
                              |_________|          0.17

    Parameters
    --------------
    Cov_result: dictionary
        Keys are virus IDs, and values are dictionaries with keys 'trim_ref_nuc', 'trim_map_nuc',
        'trim_map_loc' and values are lists of reference nucleotides, mapped nucleotides, mapped positions in reference genome.
        example: {‘NC.006432’: {‘trim_map_loc’: [3517, 3518, 3519, 3522, 3523], ‘trim_ref_nuc’: [‘A’, ‘A’, ‘C’, ‘T’, ‘A’],
         ‘trim_map_nuc’: [‘AAAA’, ‘AAT’, ‘CC’, ‘TTGTT’, ‘AAAA’]},
        ‘NC.005269’: {‘trim_map_loc’: [201, 202], ‘trim_ref_nuc’: [‘C’, ‘G’], ‘trim_map_nuc’: [‘C’, ’GG’]}}

    incoherenceThreshold: double.
        Defines maximum incoherence score allowed in order for region/positions to be kept, and
        not filtered out.

    smooth: int.
        Should be 1 or 0. 1 will make the functions smooth the incoherence score using sliding windows. 0 will make
        the function filter without any filtering.

    binsize: int, optional.
        Size of sliding window for smoothing

    jumpSize: int, optional.
        Specifies the maximum distance allowed between two mapped positions in order for the sliding window to
        keep sliding across those positions.

    Returns
    ---------------
    Dictionary with the virus ID as keys. The values are another dictionary, containing the keys
    ‘positions’, ‘mapped_nucs’, ‘ref_nuc’ and ‘nr_removed_sites’. The values of these dictionaries are lists with
    mapping locations,  mapped nucleotides, reference nucleotides and number of sites that were removed by the filter.
    """
    import numpy as np

    result = {}

    for x in Cov_result:
        Fractions = []
        Positions = []
        reference = Cov_result[x]['trim_ref_nuc']
        mapped = Cov_result[x]['trim_map_nuc']
        for m in range(len(Cov_result[x]['trim_map_nuc'])):  # Go through all virus positions present in the sample
            # (all rows in the table belonging to this virus)
            maps = Cov_result[x]['trim_map_nuc'][m]  # mapping, example  'TT'
            position = Cov_result[x]['trim_map_loc'][m]  # position, example '3517'
            a = 0
            t = 0
            c = 0
            g = 0
            for j in range(len(maps)):  # go through the mapped nucleotides

                if maps[j] == 'A':
                    a = a + 1
                elif maps[j] == 'T':
                    t = t + 1
                elif maps[j] == 'C':
                    c = c + 1
                elif maps[j] == 'G':
                    g = g + 1

            DominantNuc = max(a, t, c, g)  # Find the number of the dominant nucleotide
            DisagreeingNucs = sum([a, t, c, g]) - DominantNuc  # Find how many mapped nucs are not the dominant one
            Fraction = float(DisagreeingNucs) / float(
                sum([a, t, c, g]))  # Calculate the fraction of nucs that are not dominant
            # print [DominantNuc, DisagreeingNucs]
            Fractions.append(Fraction)
            Positions.append(position)

        FilteredFractions = []
        FilteredPositions = []
        FilteredMap = []
        FilteredRef= []
        deleteSites = []

        if smooth == 1:  # if we want to smooth before filtering
            gaps = [-1]  # save positions before a jump, ie the index before a new section starts
            pos = 0
            numlengths = []  # keep track of bin sizes

            # FIND GAPS IN MAPPING
            while pos < len(Fractions) - 1:
                if int(Positions[pos + 1]) > int(
                        Positions[pos]) + jumpSize:  # if there is a jump of at least jumpSize positions
                    gaps.append(pos)  # save the index to the positions before the jump
                pos = pos + 1

            regions = {}  # dictionary to store the sections of mapped genome
            bins = []
            newPos = []
            if len(gaps) > 1:  # check if we have jumps in the mapping
                for t in range(len(gaps)):  # Go through all the jumps
                    if gaps[t] != gaps[-1]:  # if we are not at the last jump
                        regions[t] = Positions[int(gaps[t]) + 1:int(gaps[t + 1]) + 1]  # add mapped region to dict
                    else:  # if we are at the last jump
                        regions[t] = Positions[int(gaps[t]) + 1:]  # add last mapped region to dict

                # SMOOTH
                pointer = 0  # where we are in the Fractions vector
                for key in regions:  # go through all the regions
                    nums = []  # to store values for the bin
                    if len(regions.get(key)) < binsize:  # is the region smaller than the bin?
                        for p in range(len(regions.get(key))):
                            nums.append(Fractions[pointer + p])
                        bins.append(np.mean(nums))  # take mean of the values to create the bin value
                        newPos.append(pointer)  # store the first position in the bin
                        numlengths.append(len(nums))
                        if pointer + len(nums) >= len(Positions):  # if we are at the end of the sequence
                            break
                        pointer = gaps[key + 1] + 1  # move the pointer to the start of the next region
                    else:
                        for j in range(len(regions.get(key))):  # go through positions in the region
                            nums = []
                            if j + binsize <= len(regions.get(key)):  # If we can fit another whole bin
                                for p in range(binsize):
                                    nums.append(Fractions[pointer + p])  # add values to the bin
                                pointer = pointer + 1  # move the pointer
                            else:
                                if pointer + binsize > len(Positions):  # if we are at the end of the sequence
                                    break
                                pointer = gaps[key + 1] + 1  # move the pointer to the start of the next region
                                break

                            bins.append(np.mean(nums))  # take mean of the values to create the bin value
                            newPos.append(pointer)  # store the first position in the bin
                            numlengths.append(len(nums))

            else:  # if we don't have any jumps in the mapping
                pointer = 0
                for t in range(len(Positions)):
                    nums = []
                    if pointer < len(Positions) - binsize:  # if we can fit another bin
                        for p in range(binsize):
                            nums.append(Fractions[pointer + p])
                    else:
                        break
                    bins.append(np.mean(nums))
                    newPos.append(pointer)
                    numlengths.append(len(nums))
                    pointer = pointer + 1

            for b in range(len(bins)):  # go through all the bins
                if bins[b] > incoherenceThreshold:  # Check if bin passes incoherence test
                    for d in range(numlengths[b]):
                        deleteSites.append(int(newPos[b]) + d)  # save the indices we want to delete

            deleteSites = list(set(deleteSites))  # get rid of duplicate indices
            for y in range(len(Positions)):
                if y not in deleteSites:
                    FilteredPositions.append(Positions[y])
                    FilteredFractions.append(Fractions[y])
                    FilteredRef.append(reference[y])
                    FilteredMap.append(mapped[y])
        else:
        # IF WE DONT WANT TO SMOOTH, USE THIS PART FOR FILTERING
            for b in range(len(Positions)):
                if Fractions[b] < incoherenceThreshold:
                    FilteredFractions.append(Fractions[b])
                    FilteredPositions.append(Positions[b])
                    FilteredRef.append(reference[b])
                    FilteredMap.append(mapped[b])
        result[x] = {}
        result[x]['positions'] = FilteredPositions
        result[x]['mapped_nucs'] = FilteredMap
        result[x]['ref_nuc'] = FilteredRef
        result[x]['nr_removed_sites'] = len(deleteSites)

    return result