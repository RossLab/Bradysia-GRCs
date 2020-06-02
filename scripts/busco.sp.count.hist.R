#setwd("/Users/christina/projects/Sciara-L-chromosome/")

#1. making a histogram of number of species per Busco gene type:
multicp <- read.table('data/sp.count.dup.busco.txt', stringsAsFactors = F, header = F, col.names = c("ID", "sp_count"))
singlecp <- read.table('data/sp.count.sc.busco.txt', stringsAsFactors = F, header = F, col.names = c("ID", "sp_count"))


hist(singlecp$sp_count, breaks = 10, col="green")
hist(multicp$sp_count, breaks = 10, col="blue", add=T)

#there are 16 sp. total so taking genes with 75% of sp. being present would mean any gene with 12 or more sp.
#let's see how many genes that would exclude



exclude<-singlecp[singlecp$sp_count<12,]
#57 genes
exclude2<-multicp[multicp$sp_count<12,]
#22 genes

exclude3<-singlecp[singlecp$sp_count<12.8,]
#100 genes
exclude4<-multicp[multicp$sp_count<12.8,]
#45 genes


##1345 genes 
#so that means we would exclude 79 total (5.8%) I think that's fine (for 75% of sp. present)
#80% sp. present 145 genes (10.8%)
#I think I'm going to include only genes that have more than 80%of the species
all.sp<-rbind(multicp,singlecp)
sp.list<-all.sp[all.sp$sp_count>12.8,]
sp.ID<-sp.list$ID
write.table(sp.ID, 'data/busco.gene.80.txt', quote = F, sep = "\t", row.names = F, col.names = F)

