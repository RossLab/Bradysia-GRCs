#!/usr/bin/env python3

from collections import defaultdict
import argparse
from sys import stdout
import uuid

parser = argparse.ArgumentParser(description='Reciprocal BLAST: from blast std output 6 generates two tables of orthologs as two with reciprocal hits fulfilling specified criteria')
parser.add_argument('blastfile', help='blast file (formated as 6)')
parser.add_argument('output_pattern', help='pattern for generating output')
parser.add_argument('-s', '-similarity', type = float, help='sequence similarity to accept the paralogy relation', default = 50)
parser.add_argument('-gene_coverage', type = float, help='minimal fraction of gene to be covered to accept a link of two genes', default = 0.6)

args = parser.parse_args()
args.blastfile = 'data/genome_wide_paralogy/genes_all_vs_all.blast'
args.output_pattern = 'test'
args.s = 60
args.gene_coverage = 0.6

class blast_std6_qlen_slen:
    def __init__(self, blast_line):
        blast_list = blast_line.rstrip('\n').split('\t')
        self.query = blast_list[0]
        self.subject = blast_list[1]
        self.pident = float(blast_list[2]) # percentage of identical matches
        self.len = int(blast_list[3]) # alignment length (sequence overlap)
        self.mismatch = int(blast_list[4]) # number of mismatches
        self.gapopen = int(blast_list[5]) # number of gap openings
        self.qstart = int(blast_list[6])
        self.qend = int(blast_list[7])
        self.sstart = int(blast_list[8])
        self.send = int(blast_list[9])
        self.evalue = float(blast_list[10]) # expect value
        self.bitscore = float(blast_list[11]) # bit score
        self.qlen = int(blast_list[12]) # length of query
        self.slen = int(blast_list[13]) # length of subject

    def __lt__(self, other):
        return self.qstart < other.qstart


blast_input_dict = defaultdict(list)

with open(args.blastfile, 'r') as blastf:
    # this loads all the entries where the query (entry.split('\t')[0]) and the subject (entry.split('\t')[1]) are not the same
    for entry in blastf.readlines():
        aln = blast_std6_qlen_slen(entry)
        if aln.subject != aln.query:
            blast_input_dict[aln.query].append(aln)

blast_input = list()

# curate alignments
# build dictionary of relations (every quary points to the list of subjects it's aligned to)
# query = 'jg27403'

for query in blast_input_dict.keys():
    # ALL ADJUSTMENTS of the contition for gene <-> gene orthology
    # SHOULD BE ADDED HERE
    aln_list = defaultdict(list)
    for aln in blast_input_dict[query]:
        aln_list[aln.subject].append(aln)

    for subject in aln_list.keys():
        list_of_all_q2s_alignments = aln_list[subject]
        list_of_all_q2s_alignments.sort()
        final_q2s_aln = list_of_all_q2s_alignments[0]
        for other_aln in list_of_all_q2s_alignments[1:]:
            start = min(final_q2s_aln.qstart, other_aln.qstart)
            end = max(final_q2s_aln.qend, other_aln.qend)
            # new identity is a weighted mean of ideanties of the original alginment
            final_q2s_aln.pident = ((final_q2s_aln.pident * final_q2s_aln.len) + (other_aln.pident * other_aln.len)) / (final_q2s_aln.len + other_aln.len)
            if final_q2s_aln.qend > other_aln.qstart:
                final_q2s_aln.len = end - start
            else:
                final_q2s_aln.len = final_q2s_aln.len + other_aln.len
            final_q2s_aln.qstart = start
            final_q2s_aln.qend = end
        blast_input.append(final_q2s_aln)


query2subject = defaultdict(list)

for aln in blast_input:
    # ALL ADJUSTMENTS of the contition for gene <-> gene orthology
    # SHOULD BE ADDED HERE
    if aln.pident > args.s and (aln.len / aln.qlen) > args.gene_coverage:
        query2subject[aln.query].append(aln.subject)

# query2subject['jg21908']

list_of_genes = list(query2subject.keys())
gene2og = defaultdict(str)
og2gene = defaultdict(list)
og_pairs = set()

# print('everything is ready')

for gene in list_of_genes:
    for ortholog in query2subject[gene]:
        # proceed only if the blast is reciprocal
        if not gene in query2subject[ortholog]:
            continue
        og_pairs.add(tuple(sorted((gene, ortholog))))
        # both the gene and the ortholog have OG assigned already
        if gene2og[gene] and gene2og[ortholog]:
            # do the have the same one?
            if gene2og[gene] == gene2og[ortholog]:
                continue
            # no? then merge them, geping always the gene_og, but does not really matter
            gene_og = gene2og[gene]
            for orthologs_to_merge in og2gene[gene2og[ortholog]]:
                gene2og[orthologs_to_merge] = gene_og
                og2gene[gene_og].append(orthologs_to_merge)
            continue
        # ortholog has a OG -> add the gene there as well
        if gene2og[ortholog] and not gene2og[gene]:
            gene2og[gene] = gene2og[ortholog]
            og2gene[og_id].append(gene)
            continue
        # gene has a OG -> add the ortholog there as well
        if gene2og[gene] and not gene2og[ortholog]:
            gene2og[ortholog] = gene2og[gene]
            og2gene[og_id].append(ortholog)
            continue
        # none has, make a new OG
        og_id = uuid.uuid4()
        while og_id in gene2og.values():
            og_id = uuid.uuid4()
        gene2og[ortholog] = og_id
        gene2og[gene] = og_id
        og2gene[og_id].append(ortholog)
        og2gene[og_id].append(gene)

og_hast2og_name = dict()
with open(args.output_pattern + "_OGs.tsv", 'w') as og_file:
    og_number = 1
    for og in og2gene.keys():
        og_name = 'og_' + str(og_number)
        og_hast2og_name[og] = og_name
        og_file.write("\t".join([og_name] + og2gene[og]) + '\n')
        og_number += 1

with open(args.output_pattern + "_OG_pairs.tsv", 'w') as og_pairs_file:
    og_pairs_file.write("\t".join(['og', 'gene1', 'gene2', 'identity', 'aln_len'])  + '\n')
    for aln in blast_input:
        if (aln.query, aln.subject) in og_pairs:
            og_pairs_file.write("\t".join([og_hast2og_name[gene2og[aln.query]], aln.query, aln.subject, str(aln.pident), str(aln.len),  str(aln.qlen), str(aln.slen)])  + '\n')

# for gene1, gene2 in og_pairs:
#    print("\t".join([og_hast2og_name[gene2og[gene1]], gene1, gene2]))