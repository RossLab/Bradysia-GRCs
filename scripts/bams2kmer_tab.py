# arguments L_bam, X_bam, A_bam

L_bamfile = 'data/L-27mer_mapped_to_sample.bam'
X_bamfile = 'data/X-27mer_mapped_to_sample.bam'
A_bamfile = 'data/A-27mer_mapped_to_sample.bam'

import pysam
from collections import defaultdict

for bamfile in [L_bamfile, X_bamfile, A_bamfile]:
	with pysam.AlignmentFile(bamfile, "rb") as bam:
		if not bam.check_index():
			pysam.index(bamfile)


A_bam = pysam.AlignmentFile(A_bamfile, "rb")
if not A_bam.check_index():
	pysam.index(A_bamfile)
	A_bam.close()
	A_bam = pysam.AlignmentFile(A_bamfile, "rb")

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

read_table = defaultdict(mapped_kmers)

with pysam.AlignmentFile(L_bamfile, "rb") as L_bam:
	for read in L_bam.fetch():
		read_table[read.reference_name].addL()

with pysam.AlignmentFile(A_bamfile, "rb") as A_bam:
	for read in A_bam.fetch():
		read_table[read.reference_name].addA()

with pysam.AlignmentFile(X_bamfile, "rb") as X_bam:
	for read in X_bam.fetch():
		read_table[read.reference_name].addX()

read_table_file = "data/read_mpped_kmer.tsv"
with open(read_table_file, 'w') as outtab:
	outtab.write('len\tL\tX\tA\n')
	for read in read_table.keys():
		read_mapped_kmers = read_table[read]
		start, end = read.split("/")[2].split("_")
		read_length = int(end) - int(start)
		outtab.write('{}\t{}\t{}\t{}\n'.format(read_length, read_mapped_kmers.L, read_mapped_kmers.X, read_mapped_kmers.A))