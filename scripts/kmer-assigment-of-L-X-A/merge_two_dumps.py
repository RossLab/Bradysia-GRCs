#!/usr/bin/env python3

from sys import stderr
import sys

dump_1_file = sys.argv[1]
dump_2_file = sys.argv[2]

# test for existence of the two files??

stderr.write('Merging two dump files\n')
stderr.write('Dump file 1:' + dump_1_file + '\n')
stderr.write('Dump file 2:' + dump_2_file + '\n')
stderr.write('Output: streamed to stdout\n')

i = 0
with open(dump_1_file) as d1, open(dump_2_file) as d2:
	d1_kmer = d1.readline().rstrip().split('\t')
	d2_kmer = d2.readline().rstrip().split('\t')
	while d1_kmer[0] and d2_kmer[0]:
		if d1_kmer[0] < d2_kmer[0]:
			# print(d1_kmer[0] + " < " + d2_kmer[0])
			print(d1_kmer[0] + "\t" + d1_kmer[1] + "\t0")
			d1_kmer = d1.readline().rstrip().split('\t')
		elif d1_kmer[0] > d2_kmer[0]:
			# print(d1_kmer[0] + " > " + d2_kmer[0])
			print(d2_kmer[0] + "\t0\t" + d2_kmer[1])
			d2_kmer = d2.readline().rstrip().split('\t')
		elif d1_kmer[0] == d2_kmer[0]:
			# print(d1_kmer[0] + " == " + d2_kmer[0])
			print(d1_kmer[0] + "\t" + d1_kmer[1] + "\t" + d2_kmer[1])
			d1_kmer = d1.readline().rstrip().split('\t')
			d2_kmer = d2.readline().rstrip().split('\t')
	# one one of the dumps is running out of kmers dump all from the other
	while d1_kmer[0]:
		print(d1_kmer[0] + "\t" + d1_kmer[1] + "\t0")
		d1_kmer = d1.readline().rstrip().split('\t')
	while d2_kmer[0]:
		print(d2_kmer[0] + "\t0\t" + d2_kmer[1])
		d2_kmer = d2.readline().rstrip().split('\t')
