blast_out <- read.delim('data/genome_wide_paralogy/Mdes_homologue_k05819_2_Scop_genome.blast', header=F, stringsAsFactors = F, col.names =c("query","subject","identity","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"))
scf_asn <- read.delim("data/scaffold_assignment_tab_full.tsv", header=T, stringsAsFactors = F)

row.names(scf_asn) <- scf_asn$scf

scf_asn[unique(blast_out$subject), ]

for (scf in unique(blast_out$subject)){
    one_scf <- blast_out[blast_out$subject == scf, ]
    # here I wanted to check lengths of alligments
}
