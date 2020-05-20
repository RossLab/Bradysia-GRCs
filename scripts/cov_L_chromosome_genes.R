#setwd("/Users/christina/projects/Sciara-L-chromosome/")

#1. read in gene coverage file and contig assignment file and make list of L contig ID's
scf_asn <- read.table("data/scaffold_assignment_tab_full.tsv", header = T, sep = '\t', stringsAsFactors = F)
cov_genes <- read.table("data/gene.cov.braker.annotation.tsv", header = F, sep = '\t', stringsAsFactors = F,col.names = c('scf', 'annotation.source', 'feature', 'start', 'end', 'braker1', 'braker2', 'braker3', 'braker4', 'mean_cov'))
blast_paralogs <- read.table("data/blast_paralog_table.tsv", header = T, sep = '\t', stringsAsFactors = F)


L_contig<- scf_asn[scf_asn$assignments=="L",]
L_id<-L_contig[c(2)]
#length=[1] 10628

#now I'm going to make a table with the coverage values of just the L contig genes
cov_L <- merge(cov_genes, L_id)
#length=15812
#max(cov_L$mean_cov)=34141.8

#some plots of basically the same thing
hist(cov_L$mean_cov, xlim=c(0,80),breaks=100000)


Ldens <- density(cov_L$mean_cov[cov_L$mean_cov < 80])
second_deriv <- diff(sign(diff(Ldens$y)))
plot(Ldens, xlim=c(0,80))
Ldens$x[which(second_deriv == -2) + 1]
lines(c(24.64183, 24.64183), c(0, 1), lty = 3, col = 'red')
lines(c(30.29147, 30.29147), c(0, 1), lty = 3, col = 'red')

#Now we would like to see whether paralogs are in different peaks in the histogram



paralogsLL<-blast_paralogs[blast_paralogs$subject_asn=='L',]
paralogsLL<-paralogsLL[paralogsLL$query_asn=='L',]
paralogsLL <- paralogsLL[!is.na(paralogsLL$query_asn),]

head(paralogsLL, n=30)

geneid<-cov_genes$braker4
#write.table(geneid, 'data/geneid', quote = F, sep = "\t", row.names = F)


#cov_genes$geneid_braker<-substr(cov_genes$braker4,4,10)


#gene_cov_tab<-cov_genes[c(10,11)]
blast_para<-paralogsLL[c(1,2)]

#cut -f 2 -d "=" geneid.tsv | cut -f 1 -d ";" > geneid2.tsv
gene_braker <- read.table("data/geneid2.tsv", header = T, sep = '\t', stringsAsFactors = F)
names(gene_braker)[names(gene_braker) == "x"] <- "braker_gene"
gene_cov_tab<-cbind(cov_genes,gene_braker)
gene_cov_tab<- gene_cov_tab[c(10,12)]


cov_tab <- merge(blast_para,gene_cov_tab,by.x="query",by.y = "braker_gene")
names(cov_tab)[names(cov_tab) == "mean_cov"] <- "mean_cov_query"
#length(cov_tab$query)
#length(gene_cov_tab$mean_cov)
#length(blast_para$query)

cov_tab2 <- merge(cov_tab,gene_cov_tab,by.x="subject",by.y = "braker_gene")
names(cov_tab2)[names(cov_tab2) == "mean_cov"] <- "mean_cov_subject"
#length(cov_tab2$subject)
cov_tab2$covdif<-cov_tab2$mean_cov_subject/cov_tab2$mean_cov_query

smaller_ones <- apply(cov_tab2[,c(3,4)], 1, min)
bigger_ones <- apply(cov_tab2[,c(3,4)], 1, max)

#and plot the two histograms

hist(bigger_ones, col = 'purple', breaks = 50000, xlim=c(0,100))
hist(smaller_ones, col = 'yellow', breaks=2500, xlim=c(0,100),add = T)
#hist(cov_tab2$covdif, xlim=c(0,5), breaks=100000)
hist(smaller_ones/(smaller_ones+bigger_ones), breaks=2500)
final_tab<-data.frame(smaller_ones,bigger_ones)

write.table(final_tab, 'data/mean_cov_Lparalogs.tsv', quote = F, sep = "\t", row.names = F)
