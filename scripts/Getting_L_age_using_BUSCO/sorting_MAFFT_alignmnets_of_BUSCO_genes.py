#!/usr/bin/env python3

from collections import defaultdict
import argparse
from sys import stdout
from sys import stderr
from Bio import SeqIO, Seq, SeqRecord

parser = argparse.ArgumentParser(description='Input the MAFFT alignment in fasta format')
parser.add_argument('fastafile', help='The input fasta file (all the headers are expected in the ">Genus_species_whatever" format)')
# parser.add_argument('output', default = stdout, help='pattern for generating output')

args = parser.parse_args()

odered_species = ["sylvicola_fuscatus", "penthetria_funebris", "bolitophila_cinerea", "bolitophila_hybrida", "catotricha_subobsoleta", "diadocidia_ferruginosa", "gnoriste_bilineata", "lestremia_cinerea", "macrocera_vittata", "mayetiola_destructor", "phytosciara_flavipes", "platyura_marginata", "sciara_coprophila", "symmerus_nobilis", "trichosia_splendens", "porricondyla_nigripennis"]
outgroups = ["sylvicola_fuscatus", "penthetria_funebris"]
sp2seq = defaultdict(list)

for seq_record in SeqIO.parse(args.fastafile, "fasta"):
    species_name = "_".join(seq_record.name.split('_')[0:2])
    sp2seq[species_name].append(seq_record)

outgroup = False
for sp in odered_species:
    sequences = sp2seq[sp]
    if sp == "sylvicola_fuscatus":
        if len(sequences) > 1:
            #report problem
            stderr.write(args.fastafile + ': More than one sylvicola fuscatus sequence (skipping this record).\n')
            continue
        if len(sequences) == 1:
            outgroup = True
    if sp == "penthetria_funebris" and len(sequences) == 1 and not outgroup:
        outgroup = True
    if outgroup:
        for seq in sequences:
            stdout.write(">" + seq.name + "\n")
            stdout.write(str(seq.seq) + "\n")

# with open(args.blastfile, 'r') as blastf:
#     # this loads all the entries where the query (entry.split('\t')[0]) and the subject (entry.split('\t')[1]) are not the same
#     blast_input = [entry.split('\t') for entry in blastf.readlines() if entry.split('\t')[0] != entry.split('\t')[1]]
