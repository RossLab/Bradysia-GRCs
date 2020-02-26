#!/usr/bin/env python3

from collections import defaultdict
from sys import stderr
import sys

merged_dump_file = sys.argv[1]
hist_file = sys.argv[2]

hist_tab = defaultdict(int)
max_cov = 0

with open(merged_dump_file, 'r') as dump:
	for line in dump:
		kmer = line.rstrip().split('\t')
		head = int(kmer[1])
		testes = int(kmer[2])
		if head < 5 and 10 < testes:
			hist_tab[testes] += 1
			if testes > max_cov:
				max_cov = testes

with open(hist_file, 'w') as hist:
	for cov in range(max_cov):
		hist.write(str(cov) + "\t" + str(hist_tab[cov]) + "\n")
