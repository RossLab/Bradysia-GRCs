from collections import defaultdict
from sys import stdout
import re

gene2cov = defaultdict(int)

with open('data/gene.cov.braker.annotation.tsv') as gene_coverage_file:
    for line in gene_coverage_file:
        gene_tab = line.rstrip('\n').split('\t')
        gene_name = re.match(r"ID=(.*);", gene_tab[8]).group(1)
        gene2cov[gene_name] = float(gene_tab[9])

scf2asn = dict()
with open("data/scaffold_assignment_tab_full.tsv") as scf2asn_file:
    for line in scf2asn_file:
        scf_tab = line.rstrip('\n').split('\t')
        scf = '_'.join(scf_tab[0].split('_')[0:2])
        scf2asn[scf] = scf_tab[13]

gene2asn = dict()
gene2scf = dict()
with open("data/genome/gene.scaffold.map.tsv") as gene2scf_file:
    for line in gene2scf_file:
        gene, scf = line.rstrip('\n').split('\t')
        gene2scf[gene] = scf
        gene2asn[gene] = scf2asn[scf]

with open('data/genome_wide_paralogy/anchored_filtered/Scop_anch_prot.collinearity') as colinearity_file:
    block = 0
    anch1 = ''
    anch2 = ''
    stdout.write('gene\tchr\torig_scf\tanchored_scf\tblock\torder_in_block\torig_cov\tparalog_id\tparalog_asn\n')
    for line in colinearity_file:
        if line.startswith("## Alignment"):
            block += 1
            gene_in_block = 1
            anch1, anch2 = line.split(' ')[6].split('&')
        if line.startswith("#"):
            continue

        gene1, gene2 = line.split('\t')[1:3]
        stdout.write(("{}\t" * 8 + "{}\n").format(gene1, gene2asn[gene1], gene2scf[gene1], anch1, str(block) + "a", gene_in_block, gene2cov[gene1], gene2, gene2asn[gene2]))
        stdout.write(("{}\t" * 8 + "{}\n").format(gene2, gene2asn[gene2], gene2scf[gene2], anch2, str(block) + "b", gene_in_block, gene2cov[gene1], gene1, gene2asn[gene1]))

        gene_in_block += 1
