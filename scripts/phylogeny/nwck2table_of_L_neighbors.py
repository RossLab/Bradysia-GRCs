# input_newick = "data/GRC_phylogenies/15_genetrees_L/124433at50557.treefile"
# input_nexfile = "data/GRC_phylogenies/15_genetrees_L/67606at50557.splits.nex"

import os
from collections import defaultdict
import numpy as np
from Bio import Phylo

l_string = 'L-sciara_coprophila'
a_string = 'A-sciara_coprophila'
sciaridae = set(['A-sciara_coprophila', 'phytosciara_flavipes', 'trichosia_splendens'])
cecidomyiidae = set(['mayetiola_destructor', 'porricondyla_nigripennis', 'catotricha_subobsoleta', 'lestremia_cinerea'])

def tip2sp_name(tip):
    return("_".join(tip.name.split('_')[:2]).rstrip("'"))

def path2node(path):
    if len(path) == 1:
        return(path[0])

    for clade in reversed(path[:-1]):
        terminals = clade.get_terminals()
        # if all tips are L tips, skip
        if all([tip2sp_name(t) == 'L-sciara_coprophila' for t in terminals]):
            continue
        # if the bootstrap confidence is smaller than...
        try:
            if float(clade.name.split('/')[0]) < 60:
                continue
        except AttributeError:
            continue
        break
    return(clade)

def indices2assignment(clades):
    member_sciaridae = False
    member_cecidomyiidae = False
    member_others = False
    for clade in clades:
        sp = tip2sp_name(clade)
        if sp == l_string:
            continue
        elif sp in sciaridae:
            member_sciaridae = True
            continue
        elif sp in cecidomyiidae:
            member_cecidomyiidae = True
            continue
        else:
            member_others = True
    if member_sciaridae and not member_cecidomyiidae and not member_others:
        return "sciaridae"
    if member_cecidomyiidae and not member_sciaridae and not member_others:
        return "cecidomyiidae"
    return "other"

def tree2assigments(input_newick):
    tree = Phylo.read(input_newick, "newick")

    sp2node_name = defaultdict(list)
    for tip in tree.get_terminals():
        sp = "_".join(tip.name.split('_')[:2]).rstrip("'")
        sp2node_name[sp].append(tip)

    assignments = []
    branch_lengths = []
    for l_clade in sp2node_name[l_string]:
        last_node = path2node(tree.get_path(l_clade))
        monophy_terminal = last_node.get_terminals()
        l_asn = indices2assignment(monophy_terminal)
        assignments.append(l_asn)
        branch_lengths.append(l_clade.branch_length)

    while len(assignments) < 2:
        assignments.append('NA')
        branch_lengths.append(-1)

    if len(sp2node_name[a_string]) == 1:
        a_tip = sp2node_name[a_string][0]
        last_node = path2node(tree.get_path(a_tip))
        monophy_terminal = last_node.get_terminals()
        a_asn = indices2assignment(monophy_terminal)
        assignments.append(a_asn)
        branch_lengths.append(a_tip.branch_length)
    else:
        assignments.append("other")
        branch_lengths.append(-1)

    return("\t".join(assignments + [str(l) for l in branch_lengths]))


input_newick = "data/GRC_phylogenies/15b_genetrees_LL/127380at50557.treefile"

LL_files = [i for i in os.listdir('data/GRC_phylogenies/15b_genetrees_LL') if i.endswith('treefile')]
L_files = [i for i in os.listdir('data/GRC_phylogenies/15_genetrees_L') if i.endswith('treefile')]

with open('tables/L-busco-phylogenies-summary.tsv', 'w') as tab:
    tab.write('BUSCI_id\tGRC_1\tGRC_2\tCore_genome\tGRC_1_len\tGRC_2_len\tCore_len\n')

    for file in LL_files:
        gene = file.split('.')[0]
        input_newick = 'data/GRC_phylogenies/15b_genetrees_LL/' + file
        tab.write(gene + '\t' + tree2assigments(input_newick) + '\n')

    for file in L_files:
        gene = file.split('.')[0]
        input_newick = 'data/GRC_phylogenies/15_genetrees_L/' + file
        tab.write(gene + '\t' + tree2assigments(input_newick) + '\n')


