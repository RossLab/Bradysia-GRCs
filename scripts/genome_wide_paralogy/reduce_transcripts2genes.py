#!/usr/bin/env python3

from collections import defaultdict
from Bio import SeqIO
import argparse
from sys import stdout

parser = argparse.ArgumentParser(description='Given gene2transcript map extract longest transcript per gene from an fa/faa file.')
parser.add_argument('g2t_map', help='path to gene to transcript map (2 column tsv file)')
parser.add_argument('transcriptome', help='fasta or faa file (headers must match transcripts in the map)')
parser.add_argument('-filter_list', help='a list of transcripts we want to filter out', default = '')
parser.add_argument('-keep_list', help='a list of transcripts to keep, filter_list has a priority over keep list [datault: all]', default = '')

args = parser.parse_args()

filter_dict = defaultdict(bool)
keep_dict = defaultdict(bool)

if args.filter_list:
    with open(args.filter_list, 'r') as filter_list:
        for line in filter_list:
            filter_dict[line.rstrip('\n')] = True

if args.keep_list:
    with open(args.keep_list, 'r') as keep_list:
        for line in keep_list:
            keep_dict[line.rstrip('\n')] = True

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

filter_list = defaultdict()

# read transcripts
for gene in genes2transcript_seuqneces.keys():
    transcripts = genes2transcript_seuqneces[gene]
    lengest = 0
    for index, tr in enumerate(transcripts):
        if len(transcripts[lengest].seq) > len(tr.seq):
            lengest = index
        # if lengest != 0:
        #     print(lengest)
    if filter_dict[transcripts[lengest].name]:
        continue
    if not (keep_dict[transcripts[lengest].name] and args.keep_list):
        continue
    fasta_string = '>' + gene + '\n' + str(transcripts[lengest].upper().seq) + '\n'
    stdout.write(fasta_string)
