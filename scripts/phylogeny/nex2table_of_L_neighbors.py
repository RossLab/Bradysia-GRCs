#!/bin/python3

import os
from collections import defaultdict
import numpy as np

l_string = 'L-sciara_coprophila'
sciaridae = set(['A-sciara_coprophila', 'phytosciara_flavipes', 'trichosia_splendens'])
cecidomyiidae = set(['mayetiola_destructor', 'porricondyla_nigripennis', 'catotricha_subobsoleta', 'lestremia_cinerea'])


def indices2assignment(tree_node, index2sp):
    member_sciaridae = False
    member_cecidomyiidae = False
    member_others = False
    for i in tree_node:
        if index2sp[i] == l_string:
            continue
        elif index2sp[i] in sciaridae:
            member_sciaridae = True
            continue
        elif index2sp[i] in cecidomyiidae:
            member_cecidomyiidae = True
            continue
        else:
            member_others = True
    if member_sciaridae and not member_cecidomyiidae and not member_others:
        return "sciaridae"
    if member_cecidomyiidae and not member_sciaridae and not member_others:
        return "cecidomyiidae"
    return "other"

def get_L_sister(input_nexfile):

    L_sp_indexes = set()
    index2sp = dict()
    l2tree_nodes = defaultdict(list)

    with open(input_nexfile, 'r') as nfile:
        read_species = False
        read_tree = False

        for line in nfile:
            line = line.rstrip('\n')
            if line == 'TAXLABELS':
                read_species = True
                continue
            if read_species:
                if line == ';':
                    read_species = False
                    continue
                else:
                    sp_order, sp = line.split(' ')
                    sp_order = int(sp_order[:-1][1:])
                    sp = "_".join(sp[1:].split('_')[:2]).rstrip("'")
                    index2sp[sp_order] = sp
                    if sp == l_string:
                        L_sp_indexes.add(sp_order)
            if line == 'MATRIX':
                read_tree = True
                continue
            if read_tree:
                if line == ';':
                    read_tree = False
                    break
                else:
                    line = line.split('\t')
                    support = int(line[1])
                    nodes = [int(i) for i in line[2][:-1][1:].split(' ')]
                    if len(nodes) > 1 and support > 30:
                        l_present = list(set(nodes) & L_sp_indexes)
                        if len(l_present) > 0 and not len(nodes) == len(l_present):
                            for l in l_present:
                                l2tree_nodes[l].append(nodes)

    l2minimal_nodes = dict()
    for l in l2tree_nodes:
        node_sizes = [len(nodes) for nodes in l2tree_nodes[l]]
        l2minimal_nodes[l] = l2tree_nodes[l][np.argmin(node_sizes)]

    assignments = [indices2assignment(l2minimal_nodes[l], index2sp) for l in l2minimal_nodes]
    while len(assignments) < 2:
        assignments.append('NA')

    return("\t".join(assignments))

LL_files = [i for i in os.listdir('data/GRC_phylogenies/15b_genetrees_LL') if i.endswith('nex')]
L_files = [i for i in os.listdir('data/GRC_phylogenies/15_genetrees_L') if i.endswith('nex')]

with open('tables/L-busco-phylogenies-summary.tsv', 'w') as tab:
    for file in LL_files:
        gene = file.split('.')[0]
        input_nexfile = 'data/GRC_phylogenies/15b_genetrees_LL/' + file
        tab.write(gene + '\t' + get_L_sister(input_nexfile) + '\n')

    for file in L_files:
        gene = file.split('.')[0]
        input_nexfile = 'data/GRC_phylogenies/15_genetrees_L/' + file
        tab.write(gene + '\t' + get_L_sister(input_nexfile) + '\n')

