
#!/usr/bin/env python3

from Bio import SeqIO
import sys
import gzip
from collections import defaultdict

# fasta_file = 'data/100045at50557.faa'
fasta_file = sys.argv[1]

if fasta_file[-2:] == "gz":
    ffile = SeqIO.parse(gzip.open(fasta_file, "rt"), "fasta")
else :
    ffile = SeqIO.parse(fasta_file, "fasta")

busco_genes_assigned_file = 'tables/BUSCO_assigned.tsv'
scf2asn = dict()

with open(busco_genes_assigned_file, 'r') as busco_asn_file:
    busco_asn_file.readline()
    for line in busco_asn_file:
        line = line.split('\t')
        scf = line[1].split(':')[0]
        asn = line[3]
        scf2asn[scf] = asn

sp2seq = defaultdict(list)

for seq in ffile:
    sp = '_'.join(seq.name.split('_')[0:2])
    sp2seq[sp].append(seq)

for sp in sp2seq.keys():
    # if sciara_coprofila_whatever
    if sp == 'sciara_coprophila':
        for gene in sp2seq[sp]:
            scf = "_".join(gene.id.split('_')[2:4])
            if scf2asn[scf] == 'L':
                print(">" + gene.name + "_L")
            else:
                print(">" + gene.name + "_A")
            print(gene.seq)
        continue
    the_longest = 0
    the_length = 0
    for i, gene in enumerate(sp2seq[sp]):
        if len(gene) > the_length:
            the_longest = i
            the_length = len(gene)
    seq = sp2seq[sp][i]
    print(">" + seq.name)
    print(seq.seq)
