#!/usr/bin/env python3

import os
from collections import defaultdict
import numpy as np
from Bio import Phylo
import sys

l_string = 'L-sciara_coprophila'
a_string = 'A-sciara_coprophila'
na_string = 'NA-sciara_coprophila'
# 'A-sciara_coprophila'
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
        if all(['sciara_coprophila' in tip2sp_name(t) for t in terminals]):   
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
        if sp == l_string or sp == a_string or sp == na_string:
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

    tree_type = []
    assignments = []
    branch_lengths = []
    confidence = []
    for target_clade in sp2node_name[l_string] + sp2node_name[a_string] + sp2node_name[na_string]:
        if target_clade in sp2node_name[l_string]:
            tree_type.append('L')
        if target_clade in sp2node_name[a_string]:
            tree_type.append('A')
        if target_clade in sp2node_name[na_string]:
            tree_type.append('NA')
        last_node = path2node(tree.get_path(target_clade))
        monophy_terminal = last_node.get_terminals()
        asn = indices2assignment(monophy_terminal)
        assignments.append(asn)
        branch_lengths.append(target_clade.branch_length)
        try:
            node_conf = float(last_node.name.split('/')[0])
        except ValueError:
            node_conf = -1
        confidence.append(node_conf)

    tree_type_str = ','.join(tree_type)
    asn_str = ','.join(assignments)
    branch_len_str = ','.join([str(l) for l in branch_lengths])
    confidence_str = ','.join([str(c) for c in confidence])
    return("\t".join([tree_type_str, asn_str, confidence_str, branch_len_str]))


input_dir = sys.argv[1]
tree_files = [i for i in os.listdir(input_dir) if i.endswith('treefile')]

# with open('tables/L-busco-phylogenies-summary.tsv', 'w') as tab:
sys.stdout.write('BUSCO_id\ttype\tgene_tree_location\tbootstraps\tbranch_lengths\n')

for file in tree_files:
    gene = file.split('.')[0]
    input_newick = input_dir + '/' + file
    sys.stdout.write(gene + '\t' + tree2assigments(input_newick) + '\n')


