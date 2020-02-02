#!/usr/bin/env python3

from sys import stderr
import sys

merged_dump_file = sys.argv[1]
kmer_prefix = 'data/kmers_k27'

A_file = kmer_prefix + '_A.fasta'
X_file = kmer_prefix + '_X.fasta'
L_file = kmer_prefix + '_L.fasta'

A_kmers = 0
X_kmers = 0
L_kmers = 0

with open(L_file, 'w') as L, open(X_file, 'w') as X, open(A_file, 'w') as A, open(merged_dump_file, 'r') as dump:
	for line in dump:
		kmer = line.rstrip().split('\t')
		head = int(kmer[1])
		testes = int(kmer[2])
		if 125 < head and head < 175 and 80 < testes and testes < 140:
			A_kmers += 1
			kmer_fasta_record = '>A_k' + str(A_kmers) + '\n' + kmer[0] + '\n'
			A.write(kmer_fasta_record)
		if 50 < head and head < 100 and 60 < testes and testes < 100:
			X_kmers += 1
			kmer_fasta_record = '>X_k' + str(X_kmers) + '\n' + kmer[0] + '\n'
			X.write(kmer_fasta_record)
		if head < 5 and 10 < testes:
			L_kmers += 1
			kmer_fasta_record = '>L_k' + str(L_kmers) + '\n' + kmer[0] + '\n'
			L.write(kmer_fasta_record)
