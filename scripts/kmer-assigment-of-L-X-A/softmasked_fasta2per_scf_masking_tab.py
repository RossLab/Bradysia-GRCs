import sys
from Bio import SeqIO
import gzip

# fasta_file = sys.argv[1]
fasta_file = 'data/genome/scop.spades2.min1kb.trimmed.fasta.masked'

if fasta_file[-2:] == "gz":
    ffile = SeqIO.parse(gzip.open(fasta_file, "rt"), "fasta")
else :
    ffile = SeqIO.parse(fasta_file, "fasta")

#for seq_record in SeqIO.parse(sys.argv[1], "fasta"):

sys.stdout.write("name\tseqlen\tmasked\n")

for seq_record in ffile:
    scf = seq_record.name
    scf_length = len(seq_record.seq)
    softmasked_nts = sum([nt.islower() for nt in str(seq_record.seq)])
    sys.stdout.write("{}\t{}\t{}\n".format(scf, scf_length, softmasked_nts))
