def Incoherence_filter(Cov_result, incoherenceThreshold, smooth, binsize=3, jumpSize=100):
    import numpy as np
    # This function takes as input a dictionary with virus IDs as keys, and values are dictionaries were the keys are reference
    # nucleotide, mapped nucleotides and mapped nucleotides, with the corresponding values.

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
        FilteredRef = []
        FilteredMap = []

        if smooth == 1:  # if we want to smooth before filtering
            gaps = [-1]  # save positions before a jump, ie the index before a new section starts
            pos = 0
            numlengths = []  # keep track of bin sizes

            # FIND GAPS IN MAPPING
            while pos < len(Fractions) - 1:
                if int(Positions[pos + 1]) > int(Positions[pos]) + jumpSize:  # if there is a jump of at least jumpSize positions
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
                count = 0
                pointer = 0  # where we are in the Fractions vector
                for key in regions:  # go through all the regions
                    regionpointer = 0  # where we are in the region
                    while regionpointer < len(regions.get(key)):  # go through positions in the region
                        nums = []  # to store values for the bin
                        if regionpointer < len(regions.get(key)) - binsize:  # If we can fit another whole bin
                            for p in range(binsize):
                                nums.append(Fractions[pointer + p])  # add values to the bin
                        else:
                            for r in range(len(regions.get(key)) - regionpointer):
                                nums.append(Fractions[pointer + r])  # add the last values to form the last bin
                        bins.append(np.mean(nums))  # take mean of the values to create the bin value
                        newPos.append(Positions[pointer])  # store the first position in the bin
                        numlengths.append(len(nums))
                        pointer = pointer + len(nums)  # move the pointer
                        regionpointer = regionpointer + len(nums)  # move the pointer
                        count = count + 1

            else:  # if we don't have any jumps in the mapping
                pointer = 0
                count = 0
                while pointer < len(Positions):
                    nums = []
                    if pointer < len(Positions) - binsize:
                        for p in range(binsize):
                            nums.append(Fractions[pointer + p])
                    else:
                        for r in range(len(Fractions) - pointer):
                            nums.append(Fractions[pointer + r])
                    bins.append(np.mean(nums))
                    newPos.append(Positions[pointer])
                    numlengths.append(len(nums))
                    pointer = pointer + len(nums)
                    count = count + 1


            position = 0
            for b in range(len(bins)):
                if bins[b] < incoherenceThreshold:  # Save bins with low incoherence
                    for d in range(numlengths[b]):
                        FilteredFractions.append(Fractions[position + d])  #
                        FilteredPositions.append(Positions[position + d])
                        FilteredRef.append(reference[position + d])
                        FilteredMap.append(mapped[position + d])
                position = position + numlengths[b]

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

    return result