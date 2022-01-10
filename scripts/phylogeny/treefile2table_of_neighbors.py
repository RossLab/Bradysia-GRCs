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

def tip2scf(tip):
    return("_".join(tip.name.split('_')[2:4]).rstrip("'"))

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

def get_basal_sp(last_node):
    for clade in last_node.clades:
        if clade.is_terminal():
            return(clade)
    return('NA')

def indices2assignment(clades, present_sciaridae, target_clade, basal_sp):
    member_sciaridae = False
    member_cecidomyiidae = False
    member_others = False
    # these will serve to evaluate if the gene copy falls as internal branch of sciaridae or not
    outgroup_sciaridae = present_sciaridae.copy()
    for clade in clades:
        sp = tip2sp_name(clade)
        # this handles all the cases when evaluating members of sciaridae
        if sp in present_sciaridae:
            try:
                # this is to keep track is all sciaridae are show monophyly (i.e. if the gene is _o or _i)
                outgroup_sciaridae.remove(sp)
            except KeyError:
                continue
            # bradysia nodes wont be considered for phylogenetic placement (that needs to be a different member of sciaridae)
            if 'sciara_coprophila' in sp:
                continue
            else:
                member_sciaridae = True
                continue

        if sp in cecidomyiidae:
            member_cecidomyiidae = True
            continue
        else:
            member_others = True

    if member_sciaridae and not member_cecidomyiidae and not member_others:
        if outgroup_sciaridae == set() and target_clade == basal_sp:
            return "sciaridae_o"
        else:
            return "sciaridae_i"
    if member_cecidomyiidae and not member_sciaridae and not member_others:
        return "cecidomyiidae"
    return "other"

def tree2assigments(input_newick):
    tree = Phylo.read(input_newick, "newick")

    # print('parsing: ' + input_newick)
    sp2node_name = defaultdict(list)
    present_Sciaridae = set()
    possible_Sciaridae = set(['phytosciara_flavipes', 'trichosia_splendens', 'A-sciara_coprophila', 'L-sciara_coprophila', 'NA-sciara_coprophila'])
    sciaridae_tips = []
    for tip in tree.get_terminals():
        sp = "_".join(tip.name.split('_')[:2]).rstrip("'")
        sp2node_name[sp].append(tip)
        if sp in possible_Sciaridae:
            # defining which sciaridae are present in the treefile
            present_Sciaridae.add(sp)
            sciaridae_tips.append(tip)

    # print('testing monophyly')
    # if tree.is_monophyletic(sciaridae_tips):
    #     print("sciaridae are monophyletic")
    # else:
    #     print("they are not")

    tree_type = []
    assignments = []
    branch_lengths = []
    confidence = []
    scfs = []

    for target_clade in sp2node_name[l_string] + sp2node_name[a_string] + sp2node_name[na_string]:
        if target_clade in sp2node_name[l_string]:
            tree_type.append('L')
        if target_clade in sp2node_name[a_string]:
            tree_type.append('A')
        if target_clade in sp2node_name[na_string]:
            tree_type.append('NA')
        last_node = path2node(tree.get_path(target_clade))
        basal_sp = get_basal_sp(last_node)
        monophy_terminal = last_node.get_terminals()
        asn = indices2assignment(monophy_terminal, present_Sciaridae, target_clade, basal_sp)
        if asn[0:9] == 'sciaridae' and (not tree.is_monophyletic(sciaridae_tips) or len(present_Sciaridae) != 4):
            # print('changing assignment')
            asn = 'sciaridae'
        assignments.append(asn)
        branch_lengths.append(target_clade.branch_length)
        try:
            node_conf = float(last_node.name.split('/')[0])
        except ValueError:
            node_conf = -1
        confidence.append(node_conf)
        scfs.append(tip2scf(target_clade))

    tree_type_str = ','.join(tree_type)
    asn_str = ','.join(assignments)
    branch_len_str = ','.join([str(l) for l in branch_lengths])
    confidence_str = ','.join([str(c) for c in confidence])
    scf_str = ','.join(scfs)
    return("\t".join([tree_type_str, scf_str, asn_str, confidence_str, branch_len_str]))


input_dir = sys.argv[1]
tree_files = [i for i in os.listdir(input_dir) if i.endswith('treefile')]

# with open('tables/L-busco-phylogenies-summary.tsv', 'w') as tab:
sys.stdout.write('BUSCO_id\ttype\tscfs\tgene_tree_location\tbootstraps\tbranch_lengths\n')

for file in tree_files:
    gene = file.split('.')[0]
    input_newick = input_dir + '/' + file
    sys.stdout.write(gene + '\t' + tree2assigments(input_newick) + '\n')


