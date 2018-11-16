import gzip
import os
import matplotlib.pyplot as plt
import numpy as np
# This script takes an opened file with tab separated values on the form
# VirusId	Position	ReferenceNucleotide	MappedNucleotides
# and plots incoherence of mapping

incoherenceThreshold = 0.07
incoherenceThresholdSmooth = 0.03
binsize = 3
jumpSize = 100


files = []
path = 'sample_data/'  # Directory where all files are stored
for name in os.listdir(path):  # save all files to files variable
    files.append(gzip.open(path + name))

# ApeDictionary = {}

# output = 'test_output/output.tsv' #Output path
headers = ['FileName', 'VirusID', 'VirusName', 'Abundance', 'Coverage', 'ErrorRate']

virusFile = "virus_data/virus_names.tsv"
virus = {}
with open(virusFile) as f:
    for line in f:
        (id, name) = line.strip().split("\t")
        virus[id] = name
f.close()

Apes = {}

for i in range(len(files)):  # go through all files and call functions A, B and C.
    viruslist = []
    viruses = {}
    with gzip.open(files[i].name) as apfil:
        for line in apfil:
            (seqid, loc, nuc, nucmap) = line.strip().split("\t")  # read the lines and divide by column
            viruslist.append((seqid, nuc, nucmap, loc))
    apfil.close()

    for x, y, z, a in viruslist:  # go through four elements at a time
        virus = viruses.get(x, [])
        virus.append([y, z, a])
        viruses[x] = virus

    mismatches = {}
    errorRate = {}
    for x in viruses:
        if x == 'NC_021704.1':
            counter = 0
            misnuc = 0
            Fractions = []
            Positions = []
            for m in range(len(viruses.get(x))):  # Go through all virus positions present in the sample
                # (all rows in the table belonging to this virus)
                mapping = viruses.get(x)[
                    m]  # List of mapping, example ['T', 'TT', '3257'] = [refNuc, MappedNuc, position]
                maps = mapping[1]  # The mapped nucleotides, 'TT'
                A = 0
                T = 0
                C = 0
                G = 0
                for j in range(len(maps)):  # go through the mapped nucleotides
                    counter = counter + 1
                    # if maps[j] != mapping[0]: #check if mapped nucs are same as reference
                    # misnuc = misnuc +1
                    # Check how coherent the maps are

                    if maps[j] == 'A':
                        A = A + 1
                    elif maps[j] == 'T':
                        T = T + 1
                    elif maps[j] == 'C':
                        C = C + 1
                    elif maps[j] == 'G':
                        G = G + 1

                DominantNuc = max(A, T, C, G)  # Find the number of the dominant nucleotide
                DisagreeingNucs = sum([A, T, C, G]) - DominantNuc  # Find how many mapped nucs are not the dominant one
                Fraction = float(DisagreeingNucs) / float(
                    sum([A, T, C, G]))  # Calculate the fraction of nucs that are not dominant
                # print [DominantNuc, DisagreeingNucs]
                Fractions.append(Fraction)
                Positions.append(mapping[2])

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

                #SMOOTH
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


            # filter for high incoherence with smoothed data
            FilteredFractions = []
            FilteredPositions = []
            FilteredFractionsSmooth = []
            FilteredPositionsSmooth = []

            position = 0
            for b in range(len(bins)):
                if bins[b] < incoherenceThresholdSmooth:  # Save bins with low incoherence
                    for d in range(numlengths[b]):
                        FilteredFractionsSmooth.append(Fractions[position + d])  #
                        FilteredPositionsSmooth.append(Positions[position + d])
                position = position + numlengths[b]

            # IF WE DONT WANT TO SMOOTH, USE THIS PART FOR FILTERING
            for b in range(len(Positions)):
                if Fractions[b] < incoherenceThreshold:
                    FilteredFractions.append(Fractions[b])
                    FilteredPositions.append(Positions[b])

            f, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharey=True)
            ax1.plot(Positions, Fractions)
            ax2.plot(newPos, bins)
            ax3.plot(FilteredPositions, FilteredFractions)
            ax4.plot(FilteredPositionsSmooth, FilteredFractionsSmooth)
            ax1.xaxis.set_major_formatter(plt.NullFormatter())
            ax2.xaxis.set_major_formatter(plt.NullFormatter())
            ax3.xaxis.set_major_formatter(plt.NullFormatter())
            ax4.xaxis.set_major_formatter(plt.NullFormatter())

            ax1.set_xlabel('positions')
            ax2.set_xlabel('positions')
            ax3.set_xlabel('positions')
            ax4.set_xlabel('positions')
            ax1.set_ylabel('rate of incoherence')

            ax1.set_title('original')
            ax2.set_title('Smoothed')
            ax3.set_title('Filtered, kept fraction: %f' %round(float(len(FilteredPositions))/float(len(Positions)), 2))
            ax4.set_title('Filtered Smoothed, kept fraction: %f' %round(float(len(FilteredPositionsSmooth))/float(len(Positions)), 2))


            # ax3.set_title(float(len(FilteredPositions)/len(Positions)))
            # PLOT THE GAPS
            for c in gaps:
                ax1.plot([Positions[c + 1], Positions[c + 1]], [0, 0.6], color='g', linestyle='--', linewidth=1)

            f.suptitle('file = {}, binsize = {}, S_thresh = {}, thresh = {}'.format(files[i].name, binsize, incoherenceThresholdSmooth, incoherenceThreshold))
            plt.show()

# NC_029996.1
# NC_024359.1

# THE BELOW SECTION DRAWS DIFFERENT COLOURED LINES FOR DIFFERENT SIZE JUMPS IN THE MAPPED SEQUENCE
# Herpes: NC_009334.1
#            pos = 0
#            while pos < len(Fractions)-1:
#                if int(Positions[pos + 1]) > int(Positions[pos]) + 10:  # if there is a jump of at least 10 positions
#                    if int(Positions[pos + 1]) > int(Positions[pos]) + 100:  # if there is a jump of at least 100 positions
#                        if int(Positions[pos + 1]) > int(Positions[pos]) + 1000:  # if there is a jump of at least 1000 positions
#                            plt.plot([Positions[pos], Positions[pos]], [0, 0.6], color='r', linestyle='--', linewidth=1)
#                        else:
#                            plt.plot([Positions[pos], Positions[pos]], [0, 0.6], color='y', linestyle='--', linewidth=1)
#                    else:
#                        plt.plot([Positions[pos], Positions[pos]], [0, 0.6], color='g', linestyle='--', linewidth=1)

#                pos = pos + 1