## This is a small script to generate a table that has the gene assignments for each gene and plots the numbers of genes on each chromosome
#setwd("/Users/christina/projects/Sciara-L-chromosome/")
scf_asn <- read.delim("data/scaffold_assignment_tab_full.tsv", header=T, stringsAsFactors = F)
gene_asn <- read.delim("data/gene.scaffold.map.tsv", header=F, stringsAsFactors = F, col.names = c('gene', 'scf'))
row.names(scf_asn) <- scf_asn$scf
scf_asn<-scf_asn[,c(1:3,15)]

genes.to.scf<-merge(scf_asn,gene_asn, by.x="scf", by.y= "scf")

#write.table(genes.to.scf, 'data/gene.assignment.tab.tsv', quote = F, sep = "\t", row.names = F)

### now plotting summary of number of each chromosome

gene.count<-table(genes.to.scf$assignments)
gene.count<-as.data.frame(gene.count)
colnames(gene.count)[1] <- 'chromosome'
colnames(gene.count)[2] <- 'gene_number'



chromosome.type<-c("Autosomes","X", "GRC","Unassigned")
gene.num<-as.data.frame(chromosome.type)
gene.num$chromosome.type <- factor(gene.num$chromosome.type,levels = c("Autosomes", "X", "GRC", "Unassigned"))  
gene.num$num<-c(17802,4277,15812,3527)#this doensn't include NA genes

cbbPalette <- c("#009E73", "#56B4E9","#E69F00","#999999", "#F0E442","#0072B2" , "#D55E00", "#CC79A7")
ggplot(gene.num, aes(y=num, x=chromosome.type, fill=chromosome.type)) + scale_fill_manual(values=cbbPalette) +
  geom_bar(position="dodge", stat="identity")+theme_classic()+theme(axis.text=element_text(size=12))
#+ theme(legend.text=element_text(size=12))
#41418-17802-15812-4277=3527 I'll use this number as the unassigned genes