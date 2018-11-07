import gzip
from collections import defaultdict
from functionC import functionC

sample_file = "sample_data/testapa.mpile.gz"


with gzip.open(sample_file) as apfil:
	ErrorRate = functionC(apfil)


apfil.close


