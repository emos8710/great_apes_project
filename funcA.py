import gzip

file = 'sample_data/Gorilla_beringei_graueri-Serufuli.mpile.gz'
tup = []
groups = {}
lens = []
with gzip.open(file) as f:
	for line in f:
		split = line.strip().split("\t")
		tup.append((split[0], split[3]))
#print tup
for x,y in tup:
    group = groups.get(x, [])
    group.append(y)
    groups[x] = group
for x in groups:
	print(x,groups[x])
	num_cov= len(groups[x])
	print(num_cov)
#for x in groups[x]:
#	lens.append(len(x))
#    print(lens)
#num_nuc = sum(lens)
#print(num_nuc)

	#abundance = float(num_nuc / num_cov)
	#print(num_cov,num_nuc)
    #print(abundance)


#print(len(groups['AC_000007.1']))
#print(groups['AC_000007.1'].count("A") + groups['AC_000007.1'].count("G") + groups['AC_000007.1'].count("C") + groups['AC_000007.1'].count("T"))




