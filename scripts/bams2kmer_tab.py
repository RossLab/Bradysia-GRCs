#!/usr/bin/env python3

# arguments L_bam, X_bam, A_bam

from sys import stderr
import sys
import pysam
from collections import defaultdict

L_bamfile = sys.argv[1]
X_bamfile = sys.argv[2]
A_bamfile = sys.argv[3]

stderr.write('Input files:\n')
stderr.write('	L bamfile: ' + L_bamfile + '\n')
stderr.write('	X bamfile: ' + X_bamfile + '\n')
stderr.write('	A bamfile: ' + A_bamfile + '\n')

stderr.write('Checking indexes\n')
for bamfile in [L_bamfile, X_bamfile, A_bamfile]:
	with pysam.AlignmentFile(bamfile, "rb") as bam:
		if not bam.check_index():
			stderr.write('\tIndexing ' + bamfile + '\n')
			pysam.index(bamfile)
		else:
			stderr.write('\t' + bamfile + ' is already indexed\n')


class mapped_kmers(object):
    def __init__(self):
        self.L = 0
        self.X = 0
        self.A = 0
    def __repr__(self):
        return '{}\t{}\t{}'.format(self.L, self.X, self.A)
    def __str__(self):
        return '{}\t{}\t{}'.format(self.L, self.X, self.A)
    def addL(self):
        self.L += 1
    def addX(self):
        self.X += 1
    def addA(self):
        self.A += 1

stderr.write('Processing mapping files\n')
seq_table = defaultdict(mapped_kmers)

with pysam.AlignmentFile(L_bamfile, "rb") as L_bam:
	for seq in L_bam.fetch():
		seq_table[seq.reference_name].addL()
stderr.write('\tL bamfile done\n')

with pysam.AlignmentFile(A_bamfile, "rb") as A_bam:
	for seq in A_bam.fetch():
		seq_table[seq.reference_name].addA()
stderr.write('\tA bamfile done\n')

with pysam.AlignmentFile(X_bamfile, "rb") as X_bam:
	for seq in X_bam.fetch():
		seq_table[seq.reference_name].addX()
stderr.write('\tX bamfile done\n')

seq_table_file = "table_of_mapped_kmers.tsv"
stderr.write('Writing output in +' + seq_table_file + '\n')

with open(seq_table_file, 'w') as outtab:
	outtab.write('id\tL\tX\tA\n')
	for seq in seq_table.keys():
		seq_mapped_kmers = seq_table[seq]
		outtab.write('{}\t{}\t{}\t{}\n'.format(seq, seq_mapped_kmers.L, seq_mapped_kmers.X, seq_mapped_kmers.A))

stderr.write('Everything is done now\n')