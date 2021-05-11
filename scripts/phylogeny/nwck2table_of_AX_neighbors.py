#input_newick = "15d_genetrees_AXscp/93746at50557.treefile"

import os
from collections import defaultdict
import numpy as np
from Bio import Phylo

a_string = 'A-sciara_coprophila'
sciaridae = set(['phytosciara_flavipes', 'trichosia_splendens'])
cecidomyiidae = set(['mayetiola_destructor', 'porricondyla_nigripennis', 'catotricha_subobsoleta', 'lestremia_cinerea'])

def tip2sp_name(tip):
    return("_".join(tip.name.split('_')[:2]).rstrip("'"))

def path2node(path):
    if len(path) == 1:
        return(path[0])

    for clade in reversed(path[:-1]):
        terminals = clade.get_terminals()
        # if all tips are L tips, skip
        if all([tip2sp_name(t) == 'A-sciara_coprophila' for t in terminals]):
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
        if sp == a_string:
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

    a_clade = sp2node_name[a_string][0]
    last_node = path2node(tree.get_path(a_clade))
    monophy_terminal = last_node.get_terminals()
    assignment = indices2assignment(monophy_terminal)
    branch_length = a_clade.branch_length

    # #what does this do, ask kamil (didn't change)
    # if len(sp2node_name[a_string]) == 1:
    #     a_tip = sp2node_name[a_string][0]
    #     last_node = path2node(tree.get_path(a_tip))
    #     monophy_terminal = last_node.get_terminals()
    #     a_asn = indices2assignment(monophy_terminal)
    #     assignments.append(a_asn)
    #     branch_lengths.append(a_tip.branch_length)
    # else:
    #     assignment = "other"
    #     branch_length = -1

    return(assignment + "\t" + str(branch_length))


#input_newick = "15d_genetrees_AXscp/93746at50557.treefile"

#LL_files = [i for i in os.listdir('data/GRC_phylogenies/15b_genetrees_LL') if i.endswith('treefile')]
#A_files = [i for i in os.listdir('15d_genetrees_AXscp/') if i.endswith('treefile')]
A_files = [i for i in os.listdir('7_GRC_phylogenies/15d_genetrees_AXscp/93746at50557.treefile') if i.endswith('treefile')]

with open('tables/A-busco-phylogenies-summary.tsv', 'w') as tab:
    tab.write('BUSCO_id\tCore_genome\tCore_len\n')

    for file in A_files:
        gene = file.split('.')[0]
        input_newick = '7_GRC_phylogenies/15d_genetrees_AXscp/' + file
        tab.write(gene + '\t' + tree2assigments(input_newick) + '\n')
