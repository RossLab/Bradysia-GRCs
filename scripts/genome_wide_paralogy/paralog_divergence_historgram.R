
#library(tximport)
# setwd("/Users//christina//Dropbox//Sciara//assembly.outputs//all_by_all_blast_spades")
#setwd("/Users/christina/projects/Sciara-L-chromosome/")
scf_asn <- read.delim("data/scaffold_assignment_tab_full.tsv", header=T, stringsAsFactors = F)
trans_asn <- read.delim("data/gene.scaffold.map.tsv", header=F, stringsAsFactors = F, col.names = c('gene', 'scf'))
blast_out <- read.delim("data/genes_all_vs_all.blast", header=F, stringsAsFactors = F, col.names =c("query","subject","identity","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"))
blast_out <- blast_out[blast_out$subject != blast_out$query,]

row.names(scf_asn) <- scf_asn$scf
row.names(trans_asn) <- trans_asn$gene


#making a dataframe of genes that belong to each scf/assignments and generating number of genes for each assignment
genes.to.scf<-merge(scf_asn,trans_asn)
L_genes<- genes.to.scf[genes.to.scf$assignments=="L",]
X_genes<- genes.to.scf[genes.to.scf$assignments=="X",]
A_genes<- genes.to.scf[genes.to.scf$assignments=="A",]

number.genes<-c(length(L_genes$gene),length(X_genes$gene),length(A_genes$gene))
#40939 genes total
as.numeric(number.genes)
data.frame(number.genes,type.gene=c('L','X','A'))

# trans_asn[gene, 'scf']
# scf_asn[scf, 'assigments']

blast_out$subject_asn <- scf_asn[trans_asn[blast_out$subject, 'scf'], 'assignments']
blast_out$query_asn <- scf_asn[trans_asn[blast_out$query, 'scf'], 'assignments']

#output file of blast paralogs
#write.table(blast_out, 'data/blast_paralog_table.tsv', quote = F, sep = "\t", row.names = F)


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

#now we have some all-by-all reciprocal blast results. We want to add some extra info to the blast results table
#there are several blast outputs that we generated
##1. reciprocal blast hits table with each recirocal pair of genes, similarity, and length of alignment
##2. paralog groups with all reciprocal blast hits that contains all the genes that are paralogs to each other
#we did this for aa and nt sequences and have tables at different similarity thresholds for the aa blasts

#first making a table with gene ident, length, assignment, etc...(need to get gene info from braker list)

genes.info<-genes.to.scf[,c(16,15,1)]

gene.len<- read.table("data/gene.length.txt", header=F, sep=' ', stringsAsFactors = F, col.names = c('gene','nt_length'))
aa.len <- read.table("data/aa.length.txt", header=F, sep=' ', stringsAsFactors = F, col.names = c('gene','aa_length'))

gene.info <- merge(genes.info,gene.len,by.x="gene",by.y = "gene")
gene.info_aa <- merge(genes.info,aa.len,by.x="gene",by.y = "gene")
#note: the length of gene.info_aa is one less than genes.info or gene.info (41417)


#now loading some of the blast results and merging data.frames into one table
blast_nt <- read.delim("data/nt_orthology_OG_pairs.tsv", header=T, stringsAsFactors = F)
paralog_table_full<- merge(blast_nt ,gene.info,by.x="gene1",by.y = "gene")
#> length(paralog_table_full$gene1)
#[1] 9150

colnames(paralog_table_full)[6] <- 'assignment_gene1'
colnames(paralog_table_full)[7] <- 'scf_gene1'
colnames(paralog_table_full)[8] <- 'len_gene1'

paralog_table_full<- merge(paralog_table_full ,gene.info,by.x="gene2",by.y = "gene")
colnames(paralog_table_full)[9] <- 'assignment_gene2'
colnames(paralog_table_full)[10] <- 'scf_gene2'
colnames(paralog_table_full)[11] <- 'len_gene2'
#> length(paralog_table_full$gene2)
#[1] 9150         looks ok, the length is right for this table so now I want to reorder

paralog_table_full<-paralog_table_full[, c(3,2,1,4,5,6:11)]
##ok, it all looks good
