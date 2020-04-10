#!/usr/bin/env python3

from collections import defaultdict
from Bio import SeqIO
import argparse
from sys import stdout

parser = argparse.ArgumentParser(description='Given gene2transcript map extract longest transcript per gene from an fa/faa file.')
parser.add_argument('g2t_map', help='path to gene to transcript map (2 column tsv file)')
parser.add_argument('transcriptome', help='fasta or faa file (headers must match transcripts in the map)')

args = parser.parse_args()

transcripts2genes = dict()
with open(args.g2t_map, 'r') as g2t_file:
    for line in g2t_file:
        gene, trainscript = line.rstrip('\n').split('\t')
        transcripts2genes[trainscript] = gene

genes2transcript_seuqneces = defaultdict(list)
ffile = SeqIO.parse(args.transcriptome, 'fasta')
for seq_record in ffile:
    gene = transcripts2genes[seq_record.name]
    genes2transcript_seuqneces[gene].append(seq_record)

# read transcripts
for gene in genes2transcript_seuqneces.keys():
    transcripts = genes2transcript_seuqneces[gene]
    lengest = 0
    for index, tr in enumerate(transcripts):
        if len(transcripts[lengest].seq) > len(tr.seq):
            lengest = index
        # if lengest != 0:
        #     print(lengest)
    fasta_string = '>' + gene + '\n' + str(transcripts[lengest].upper().seq) + '\n'
    stdout.write(fasta_string)
