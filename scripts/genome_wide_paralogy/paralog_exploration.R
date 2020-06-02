
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

#### nt_blast results
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
#write.table(paralog_table_full, 'data/ntgene_blast_pair_summary.tsv', quote = F, sep = "\t", row.names = F)

#doing some preliminary histograms to get an idea of where we should set length cutoffs
hist(paralog_table_full$aln_len, breaks = 500, xlim=c(0,200))
#I think setting a length cutoff to 100 makes sense from histogram
paralog_table_full$percentaln<-paralog_table_full$aln_len/paralog_table_full$len_gene1
hist(paralog_table_full$percentaln, breaks=1000)

#trying to figure out how many times each og group is in table
ogtable<-table(paralog_table_full$og)
hist(ogtable, breaks=20, ylim=c(0,500))
hist(ogtable, breaks=30)
#2609 with two paralogs (out of 4588 og's, 57% have only one paralog)
ogtable2<-as.data.frame(ogtable)
count1<-ogtable2[ogtable2$Freq==1,]
paralog_table_full<-merge(paralog_table_full , ogtable2 ,by.x="og",by.y = "Var1")
#9150, 4588
paralog_table_full <- paralog_table_full[!is.na(paralog_table_full$identity),]

#Now I want to get the coverage levels of each of these genes and compare them
cov_genes <- read.table("data/gene_cov_table.tsv", header = T, sep = '\t', stringsAsFactors = F)
#col.names = c('scf', 'annotation.source', 'feature', 'start', 'end', 'braker1', 'braker2', 'braker3', 'braker4', 'mean_cov'))
#mean_cov, braker_gene
cov_genes<-cov_genes[c(10,11)]

##I'm going to put the cov values in the full table, so I can subset differently later if I want
paralog_cov_tab <- merge(paralog_table_full,cov_genes,by.x="gene1",by.y = "braker_gene")
colnames(paralog_cov_tab)[14] <- 'meancov_gene1'
paralog_cov_tab <- merge(paralog_cov_tab,cov_genes,by.x="gene2",by.y = "braker_gene")
colnames(paralog_cov_tab)[15] <- 'meancov_gene2'
#write.table(paralog_cov_tab, 'data/ntgene_blast_pair_cov.tsv', quote = F, sep = "\t", row.names = F)




paralog2_50.100<-paralog_cov_tab[paralog_cov_tab$Freq==1,]
paralog2_50.100<-paralog2_50.100[paralog2_50.100$percentaln>0.5,]
paralog2_50.100<-paralog2_50.100[paralog2_50.100$aln_len>100,]
#2609 length
#after filtering 1718 (min 50% aln percent,min100 aln length, only freq=1)

onlyLL<-paralog2_50.100[paralog2_50.100$assignment_gene1=="L",]
onlyLL<-onlyLL[onlyLL$assignment_gene2=="L",]
onlyLL <- onlyLL[!is.na(onlyLL$identity),]
#292
smaller_ones <- apply(onlyLL[,c(14,15)], 1, min)
bigger_ones <- apply(onlyLL[,c(14,15)], 1, max)

#and plot the two histograms

hist(bigger_ones, col = 'purple', breaks = 1000, xlim=c(0,100), ylim = c(0,50))
hist(smaller_ones, col = 'yellow', breaks=100, xlim=c(0,100),add = T)
#hist(cov_tab2$covdif, xlim=c(0,5), breaks=100000)

#Now I want to see how many have two L paralogs and one of something else
paralog2_50.100<-paralog_cov_tab[paralog_cov_tab$Freq==3,]
paralog2_50.100<-paralog2_50.100[paralog2_50.100$percentaln>0.5,]
paralog2_50.100<-paralog2_50.100[paralog2_50.100$aln_len>100,]
#2609 length
#after filtering 1718 (min 50% aln percent,min100 aln length, only freq=1)

onlyLL<-paralog2_50.100[paralog2_50.100$assignment_gene1=="L",]
onlyLL<-onlyLL[onlyLL$assignment_gene2=="L",]
onlyLL <- onlyLL[!is.na(onlyLL$identity),]

head (paralog_cov_tab[paralog_cov_tab$Freq==2,], n=30)

len_dif<-paralog_cov_tab$len_gene1/paralog_cov_tab$len_gene2
hist(len_dif, breaks = 1000, xlim = c(0,20))



#hmm, need to think about, I think there are some weird things going on with the paralogs
#maybe could look at most conservative aa one? not sure if that would be dif than nt.





















#### aa_blast results (70%)
blast_aa <- read.delim("data/aa_orthology_70_OG_pairs.tsv", header=T, stringsAsFactors = F)
aa_paralog_table_full<- merge(blast_aa ,gene.info_aa,by.x="gene1",by.y = "gene")
#> length(paralog_table_full$gene1)
#[1] 9150

colnames(aa_paralog_table_full)[6] <- 'assignment_gene1'
colnames(aa_paralog_table_full)[7] <- 'scf_gene1'
colnames(aa_paralog_table_full)[8] <- 'len_gene1'

aa_paralog_table_full<- merge(aa_paralog_table_full ,gene.info_aa,by.x="gene2",by.y = "gene")
colnames(aa_paralog_table_full)[9] <- 'assignment_gene2'
colnames(aa_paralog_table_full)[10] <- 'scf_gene2'
colnames(aa_paralog_table_full)[11] <- 'len_gene2'
#> length(paralog_table_full$gene2)
#[1] 9150         looks ok, the length is right for this table so now I want to reorder

aa_paralog_table_full<-aa_paralog_table_full[, c(3,2,1,4,5,6:11)]

#doing some preliminary histograms to get an idea of where we should set length cutoffs
hist(aa_paralog_table_full$aln_len, breaks = 1000, xlim=c(0,2000))
#I think setting a length cutoff to 100 makes sense from histogram
aa_paralog_table_full$percentaln<-aa_paralog_table_full$aln_len/aa_paralog_table_full$len_gene1
hist(aa_paralog_table_full$percentaln, breaks=500)

ogtable_aa<-table(aa_paralog_table_full$og)
hist(ogtable_aa, breaks=200)
hist(ogtable_aa, breaks=250, xlim=c(0,50))
#lots more (than nt) with a lot of paralogs but maybe twice the amount max with only two paralogs (3000/5679)