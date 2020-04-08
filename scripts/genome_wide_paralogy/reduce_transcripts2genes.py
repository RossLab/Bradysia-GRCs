#!/usr/bin/env python3

from collections import defaultdict
from Bio import SeqIO

transcripts2genes = dict()
gene2transcript_mapfilename = 'data/genome/transcripts2genes.map'
with open(gene2transcript_mapfilename, 'r') as g2t_file:
    for line in g2t_file:
        gene, trainscript = line.rstrip('\n').split('\t')
        transcripts2genes[trainscript] = gene

genes2transcript_seuqneces = defaultdict(list)
transcripts_fasta_filename = 'data/genome/transcripts.fasta'
ffile = SeqIO.parse(transcripts_fasta_filename, 'fasta')
for seq_record in ffile:
    gene = transcripts2genes[seq_record.name]
    genes2transcript_seuqneces[gene].append(seq_record)

# read transcripts
genes_fasta_filename = 'data/genome/genes.fasta'
with open(genes_fasta_filename, 'w') as out_file:
    for gene in genes2transcript_seuqneces.keys():
        transcripts = genes2transcript_seuqneces[gene]
        lengest = 0
        for index, tr in enumerate(transcripts):
            if len(transcripts[lengest].seq) > len(tr.seq):
                lengest = index
            # if lengest != 0:
            #     print(lengest)
        fasta_string = '>' + gene + '\n' + str(transcripts[lengest].upper().seq) + '\n'
        out_file.write(fasta_string)
