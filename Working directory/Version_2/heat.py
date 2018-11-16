import gzip
import os
import csv
import pandas as pd
import numpy as np
import seaborn as sns
from collections import defaultdict

inputFile = '../../test_output/output.tsv'

data = defaultdict(list)
with open(inputFile) as f:
	for line in f:
		(fileName, virusId, virusName, abundance, coverage, errorRate) = line.strip().split("\t")
		data[fileName].append((virusName, abundance))
f.close()
print data

