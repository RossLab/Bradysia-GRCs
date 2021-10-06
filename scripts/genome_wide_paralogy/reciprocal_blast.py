#!/usr/bin/env python3

from collections import defaultdict
import argparse
from sys import stdout
from sys import stderr

parser = argparse.ArgumentParser(description='Reciprocal BLAST: from blast std output 6 generates two tables of orthologs as two with reciprocal hits fulfilling specified criteria')
parser.add_argument('blastfile', help='blast file formated as -outfmt "6 std qlen slen"')
parser.add_argument('output_pattern', help='pattern for generating output')
parser.add_argument('-s', '-similarity', type = float, help='sequence similarity to accept the paralogy relation', default = 50)
parser.add_argument('-gene_coverage', type = float, help='minimal fraction of gene to be covered to accept a link of two genes', default = 0.6)

args = parser.parse_args()
# args.blastfile = 'data/genome_wide_paralogy/proteins_all_vs_all.blast'
# args.output_pattern = 'data/genome_wide_paralogy/aa_homology_test_85'
# args.s = 85
# args.gene_coverage = 0.6

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
# blast_input_dict['jg28688'], aims at jg5502 and jg30732

# curate alignments by
#    merging all alignments into single per a pairt of subject and query
#    building dictionary of alignments (q, s) -> aln
# query = 'jg28688' - good testing one, has 2 hits jg5502 and jg30732
blast_input = dict()
for query in blast_input_dict.keys():
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
        blast_input[(aln.query, subject)] = final_q2s_aln

# all three following should be in the disctionary
# blast_input[('jg28688', 'jg5502')]
# blast_input[('jg5502', 'jg28688')]
# and all the combinaiton with jg30732

del blast_input_dict

# now build a unidirection graph of relations given
#     similarity condition
#     align / gene length condition
query2subject = defaultdict(list)
for aln in blast_input.values():
    ################################################################
    # ALL ADJUSTMENTS of the contition for gene <-> gene orthology #
    # SHOULD BE ADDED HERE                                         #
    ################################################################
    if aln.pident > args.s and (aln.len / aln.qlen) > args.gene_coverage:
        query2subject[aln.query].append(aln.subject)

# query2subject['jg28688']
# query2subject['jg5502']
# query2subject['jg30732']

# setting up variables to

# keep track of nodes with outther braches explored
gene2processed = defaultdict(bool)
# a list for genes that will be processed per orthogroup
queue_of_genes = list()
# OG number
og = 1
# again, tested on the popular trio of orthologs
# stderr.write('You forgot to comment follwing line\n')
# queue_of_genes = ['jg28688']

og_file = open(args.output_pattern + "_OGs.tsv", 'w')
og_pairs_file = open(args.output_pattern + "_OG_pairs.tsv", 'w')
og_pairs_file.write("\t".join(['og', 'gene1', 'gene2', 'identity', 'aln_len', 'qlen', 'slen'])  + '\n')

# print('everything is ready')
# go through all the genes
for top_bot_gene_iterator in list(query2subject.keys()):
    # but skip all those that were processed already
    if gene2processed[top_bot_gene_iterator]:
        continue

    # set up the queue of genes to be processed in this potential orthogroup which starts with that one gene
    # queue_of_genes = ['jg28688']
    queue_of_genes = [top_bot_gene_iterator]
    # here I will store all the info about the orthogroup
    og2gene = list()
    og_pairs = set()
    # now go through queue_of_genes list till it's empty
    while not queue_of_genes == []:
        # stderr.write(gene + '\n')
        # pop will always pop(0) the 1st element to gene
        gene = queue_of_genes.pop(0)
        og2gene.append(gene)
        gene2processed[gene] = True
        # and here I want to process all it's orthologs
        for ortholog in query2subject[gene]:
            # proceed only if the blast is reciprocal
            if not gene in query2subject[ortholog]:
                # stderr.write(ortholog + ' is not reciprocal\n')
                continue
            # do not explore already explored branches
            if gene2processed[ortholog]:
                # stderr.write(ortholog + ' is already processed\n')
                continue
            # this gene needs to be process within this orthogroup (unless it's already queued for processing)
            if not ortholog in queue_of_genes:
                # stderr.write('Adding ' + ortholog + ' to the queue\n')
                queue_of_genes.append(ortholog)

            og_pairs.add(tuple(sorted((gene, ortholog))))

    # record the group only if it contains at least two genes
    if len(og2gene) > 1:
        og_name = 'og_' + str(og)
        og_file.write("\t".join([og_name] + og2gene) + '\n')

        for g1, g2 in og_pairs:
            aln = blast_input[g1, g2]
            aln2 = blast_input[g2, g1]
            aln.pident = round(((aln.pident * aln.len) + (aln2.pident * aln2.len)) / (aln.len + aln2.len), 3)
            aln.len = (round(aln.len + aln2.len) / 2)
            og_pairs_file.write("\t".join([og_name, aln.query, aln.subject, str(aln.pident), str(aln.len),  str(aln.qlen), str(aln.slen)])  + '\n')

        og += 1

og_file.close()
og_pairs_file.close()

# for gene1, gene2 in og_pairs:
#    print("\t".join([og_hast2og_name[gene2og[gene1]], gene1, gene2]))
