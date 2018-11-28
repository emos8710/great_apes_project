def incoherence_filter(Cov_result, incoherenceThreshold, smooth, binsize=3, jumpSize=100):
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


        FilteredFractions = list(Fractions)
        FilteredPositions = list(Positions)
        FilteredRef = list(reference)
        FilteredMap = list(mapped)

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
                pointer = 0  # where we are in the Fractions vector
                for key in regions:  # go through all the regions
                    for t in range(len(regions.get(key))):  # go through positions in the region
                        nums = []  # to store values for the bin
                        if t + binsize <= len(regions.get(key)):  # If we can fit another whole bin
                            for p in range(binsize):
                                nums.append(Fractions[pointer + p])  # add values to the bin

                        else:
                            if pointer + binsize > len(Positions):  # if we are at the end of the sequence
                                break
                            pointer = gaps[key + 1] + 1  # mov the pointer to the start of the next region
                            break
                        bins.append(np.mean(nums))  # take mean of the values to create the bin value
                        newPos.append(Positions[pointer])  # store the first position in the bin
                        pointer = pointer + 1  # move the pointer


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
                    newPos.append(Positions[pointer])
                    pointer = pointer + 1

            deleteSites=[]
            for b in range(len(bins)):  # go through all the bins

                if bins[b] > incoherenceThreshold:  # Check if bin passes incoherence test
                    deleteSites = deleteSites + range(b, b + binsize)  # save the indices we want to delete
            deleteSites = list(set(deleteSites))  # get rid of duplicate indices
            for index in sorted(deleteSites, reverse=True):
                del FilteredPositions[index]  # delete the indices that belong to bins with too high incoherence
                del FilteredFractions[index]
                del FilteredRef[index]
                del FilteredMap[index]

            print 'len Pos', len(Positions)
            print 'len filtpos', len(FilteredPositions)
            print 'len delete', len(deleteSites)

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