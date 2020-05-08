#!/usr/bin/env python3

from collections import defaultdict
import argparse
from sys import stdout
import uuid

parser = argparse.ArgumentParser(description='Reciprocal BLAST: from blast std output 6 generates two tables of orthologs as two with reciprocal hits fulfilling specified criteria')
parser.add_argument('blastfile', help='blast file (formated as 6)')
parser.add_argument('output_pattern', help='pattern for generating output')
parser.add_argument('-s', '-similarity', type = float, help='sequence similarity to accept the paralogy relation', default = 50)

args = parser.parse_args()

with open(args.blastfile, 'r') as blastf:
    # this loads all the entries where the query (entry.split('\t')[0]) and the subject (entry.split('\t')[1]) are not the same
    blast_input = [entry.split('\t') for entry in blastf.readlines() if entry.split('\t')[0] != entry.split('\t')[1]]

query2subject = defaultdict(list)

# build dictionary of relations (every quary points to the list of subjects it's aligned to)

for aln in blast_input:
    # ALL ADJUSTMENTS of the contition for gene <-> gene orthology
    # SHOULD BE ADDED HERE
    if float(aln[2]) > args.s:
        query2subject[aln[0]].append(aln[1])

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
        if (aln[0], aln[1]) in og_pairs:
            og_pairs_file.write("\t".join([og_hast2og_name[gene2og[aln[0]]], aln[0], aln[1], aln[2], aln[3]])  + '\n')

# for gene1, gene2 in og_pairs:
#    print("\t".join([og_hast2og_name[gene2og[gene1]], gene1, gene2]))