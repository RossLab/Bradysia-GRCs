from collections import defaultdict

class mapped_kmers:
    def __init__(self):
        self.L = 0
        self.X = 0
        self.A = 0
    def __repr__(self):
        return '{}\t{}\t{}'.format(self.L, self.X, self.A)
    def __str__(self):
        return '{}\t{}\t{}'.format(self.L, self.X, self.A)
    def addX(self):
        self.X += 1
    def addL(self):
        self.L += 1
    def addA(self):
        self.A += 1

read_kmers = defaultdict(mapped_kmers)

L_reads = ['READ1','READ1','READ1','READ1','READ4','READ7','READ8']
X_reads = ['READ2','READ2','READ2','READ3','READ4','READ5','READ5']
A_reads = ['READ3','READ3','READ3','READ3','READ4','READ5','READ6']

for read in L_reads:
    read_kmers[read].addL()

for read in X_reads:
    read_kmers[read].addX()

for read in A_reads:
    read_kmers[read].addA()

read_kmers