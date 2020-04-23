
#library(tximport)
# setwd("/Users//christina//Dropbox//Sciara//assembly.outputs//all_by_all_blast_spades")
scf_asn <- read.delim("data/scaffold_assignment_tab_full.tsv", header=T, stringsAsFactors = F)
trans_asn <- read.delim("data/genome/gene.scaffold.map.tsv", header=F, stringsAsFactors = F, col.names = c('gene', 'scf'))
blast_out <- read.delim("data/genome_wide_paralogy/genes_all_vs_all.blast", header=F, stringsAsFactors = F, col.names =c("query","subject","identity","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"))
blast_out <- blast_out[blast_out$subject != blast_out$query,]

row.names(scf_asn) <- scf_asn$scf
row.names(trans_asn) <- trans_asn$gene

# trans_asn[gene, 'scf']
# scf_asn[scf, 'assigments']

blast_out$subject_asn <- scf_asn[trans_asn[blast_out$subject, 'scf'], 'assignments']
blast_out$query_asn <- scf_asn[trans_asn[blast_out$query, 'scf'], 'assignments']

nrow(blast_out)
table(paste(blast_out$query_asn, blast_out$subject_asn))

#Some plots and summary tables showing results

asn_table<-table(blast_out$subject_asn, blast_out$query_asn)
identities_of_LLs <- blast_out[paste0(blast_out$subject_asn, blast_out$query_asn) == 'LL', 'identity']
identities_of_XXs <- blast_out[paste0(blast_out$subject_asn, blast_out$query_asn) == 'XX', 'identity']
identities_of_AAs <- blast_out[paste0(blast_out$subject_asn, blast_out$query_asn) == 'AA', 'identity']
identities_of_LAs <- blast_out[paste0(blast_out$subject_asn, blast_out$query_asn) %in% c('LA', 'AL'), 'identity']
identities_of_LXs <- blast_out[paste0(blast_out$subject_asn, blast_out$query_asn) %in% c('LX', 'XL'), 'identity']
identities_of_AXs <- blast_out[paste0(blast_out$subject_asn, blast_out$query_asn) %in% c('AX', 'XA'), 'identity']

pal <- c('blue', 'yellow', 'green', 'orange', 'purple', 'red')

png('figures/paralog_divergence_hist.png')
hist(blast_out$identity, breaks = 100, main = 'nt identity of genome-wide paralogs', xlab = 'Identity [ % ]')
hist(identities_of_LAs, col = 'blue', add = T, breaks = 100)
hist(identities_of_LLs, col = 'yellow', add = T, breaks = 100)
hist(identities_of_AAs, col = 'green', add = T, breaks = 100)
hist(identities_of_LXs, col = 'orange', add = T, breaks = 50)
hist(identities_of_AXs, col = 'purple', add = T, breaks = 100)
hist(identities_of_XXs, col = 'red', add = T, breaks = 100)

legend('topleft', bty = 'n', c('all', 'LA', 'LL', 'AA', 'LX', 'AX', 'XX'), col = c('black', pal), pch = 20)
dev.off()