from collections import defaultdict


filter = defaultdict(bool)

with open('data/genome_wide_paralogy/list_of_TE_candidates_to_filter.tsv') as file:
    for line in file:
        filter[line.rstrip('\n')] = True

with open('data/genome_wide_paralogy/anchored/Scop_anch_prot.gff') as file:
    for line in file:
        gene = line.rstrip('\n').split('\t')
        if not filter[gene[1]]:
            print(line.rstrip('\n'))
