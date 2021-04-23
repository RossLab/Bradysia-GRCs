# Now we have some all-by-all reciprocal blast results. We want to add some extra info to the blast results table
### there are several blast outputs that we generated
### 1. reciprocal blast hits table with each recirocal pair of genes, similarity, and length of alignment, and gene lengths
### 2. paralog groups with all reciprocal blast hits that contains all the genes that are paralogs to each other
#### we did this for aa and nt sequences and have tables at different similarity thresholds for the aa blasts

## BELOW I'm looking at just the nt blast results
#### note: the blast parameters were set so only hits with an e-value < 1e-10 are kept

### First making a table with gene ident, length, assignment, etc...(need to get gene info from braker list)

#setwd("/Users/christina/projects/Sciara-L-chromosome/")
library(ggplot2)
scf_asn <- read.delim("data/scaffold_assignment_tab_full.tsv", header=T, stringsAsFactors = F)
trans_asn <- read.delim("data/gene.scaffold.map.tsv", header=F, stringsAsFactors = F, col.names = c('gene', 'scf'))
row.names(scf_asn) <- scf_asn$scf
row.names(trans_asn) <- trans_asn$gene
genes.to.scf<-merge(scf_asn,trans_asn)

genes.info<-genes.to.scf[,c(17,16,15,1)]

### now loading some of the blast results and merging data.frames into one table

#### nt_blast results
blast_nt <- read.delim("data/nt_orthology_OG_pairs.tsv", header=T, stringsAsFactors = F, col.names = c("og", "gene1" , "gene2", "identity", "aln_length", "len_gene1", "len_gene2"))
paralog_table_full<- merge(blast_nt ,genes.info,by.x="gene1",by.y = "gene")
#paralog_table_full<-paralog_table_full[,c(1:9)]
colnames(paralog_table_full)[9] <- 'assignment_gene1'
colnames(paralog_table_full)[10] <- 'scf_gene1'
colnames(paralog_table_full)[8] <- 'cov_scf1'

paralog_table_full<- merge(paralog_table_full ,genes.info,by.x="gene2",by.y = "gene")
#paralog_table_full<-paralog_table_full[,c(1:11)]
colnames(paralog_table_full)[12] <- 'assignment_gene2'
colnames(paralog_table_full)[13] <- 'scf_gene2'
colnames(paralog_table_full)[11] <- 'cov_scf2'

paralog_table_full<-paralog_table_full[, c(3,2,1,4,5,6,8,9,10,7,11,12,13)]

##ok, it all looks good, there are 3444 paralogs sets total
#write.table(paralog_table_full, 'data/ntgene_blast_pair_summary.tsv', quote = F, sep = "\t", row.names = F)

### Now I want to add the information about gene cov and the percent of the gene the alignment covers into the table and then set some cutoff threshold.
#### At the moment I think I'll set a min length cutoff of 100nt for each gene and I'll make two tables, one for paralogs where the alignment covers 80% of the gene and one where it's less.
#### Also, I'll also separate genes into two lists based on whether the blast hit is on the same scaffold or a different scaffold 

paralog_table_full$percentaln1<-paralog_table_full$aln_len/paralog_table_full$len_gene1
paralog_table_full$percentaln2<-paralog_table_full$aln_len/paralog_table_full$len_gene2
#hist(paralog_table_full$percentaln, breaks=1000)

paralog_table_full$lendif1<-paralog_table_full$len_gene1/paralog_table_full$len_gene2
paralog_table_full$lendif2<-paralog_table_full$len_gene2/paralog_table_full$len_gene1
#len=3444

#len=3399

#trying to figure out how many times each og group is in table
ogtable<-table(paralog_table_full$og)
hist(ogtable, breaks=30)
ogtable2<-as.data.frame(ogtable)
paralog_table_full<-merge(paralog_table_full , ogtable2 ,by.x="og",by.y = "Var1")

#### Now I want to get the coverage levels of each of these genes and compare them
cov_genes <- read.table("data/gene_cov_table.tsv", header = T, sep = '\t', stringsAsFactors = F)
#col.names = c('scf', 'annotation.source', 'feature', 'start', 'end', 'braker1', 'braker2', 'braker3', 'braker4', 'mean_cov'))
#mean_cov, braker_gene
cov_genes<-cov_genes[c(10,11)]

##I'm going to put the cov values in the full table, so I can subset differently later if I want
paralog_cov_tab <- merge(paralog_table_full,cov_genes,by.x="gene1",by.y = "braker_gene")
colnames(paralog_cov_tab)[19] <- 'meancov_gene1'
paralog_cov_tab <- merge(paralog_cov_tab,cov_genes,by.x="gene2",by.y = "braker_gene")
colnames(paralog_cov_tab)[20] <- 'meancov_gene2'
colnames(paralog_cov_tab)[18] <- 'og_freq1'
#paralog_cov_tab <- paralog_cov_tab[!is.na(paralog_cov_tab$identity),]
#write.table(paralog_cov_tab, 'data/ntgene_recip_blast_cov.tsv', quote = F, sep = "\t", row.names = F)


#some histograms to get an idea of what I've got

hist(paralog_cov_tab$identity)

hist(paralog_cov_tab$percentaln1, breaks = 100)
hist(paralog_cov_tab$percentaln2, breaks = 100)

hist(sqrt(paralog_cov_tab$len_gene1), breaks = 100)
hist(sqrt(paralog_cov_tab$len_gene2), breaks = 100, add=T, col="blue")
#for some reason the hist look very different for gene 1 and 2, is this because the way blast works or the way Kamil organized the table or what
hist(paralog_cov_tab$Freq, breaks = 30)
hist(paralog_cov_tab$Freq, breaks = 30, xlim = c(0,10))


hist(paralog_cov_tab$aln_length, breaks = 100)
hist(paralog_cov_tab$aln_length, breaks = 100, xlim = c(0,2000))

#### I think I might also want to split up the lists based on frequency (which is the number of paralogs), I'm not sure it I should do 3 and less or only 1 and more than one)

table(paralog_cov_tab$og_freq1)
###1    2    3     4    5    6    7    8    9   10   11   12   14   16   17   33   38 
#2184  286  507   60   30   84   49   32    9   40   33   12   14   16   17   33   38 
#2184  143  169   15    6    14    7   4    9   4     3   1     1   1     1    1   1

#with new filtering of genes at twice the length of partner
###1    2    3    4    5    6    7    8    9   10   11   12   13   16   17   21   33 
#2176  286  501   60   30   78   42   32    9   40   33   12   13   16   17   21   33

#I think the most interesting ones are the ones with one paralog (since there are by far the most of those) and the ones with three, I'm going to work with the lists separately for now

### Below is some data exploration

#### first a histogram of types of paralogs and identities for the whole table 
identities_of_LLs <- paralog_cov_tab[paste0(paralog_cov_tab$assignment_gene1, paralog_cov_tab$assignment_gene2) == 'LL', 'identity']
identities_of_XXs <- paralog_cov_tab[paste0(paralog_cov_tab$assignment_gene1, paralog_cov_tab$assignment_gene2) == 'XX', 'identity']
identities_of_AAs <- paralog_cov_tab[paste0(paralog_cov_tab$assignment_gene1, paralog_cov_tab$assignment_gene2) == 'AA', 'identity']
identities_of_LAs <- paralog_cov_tab[paste0(paralog_cov_tab$assignment_gene1, paralog_cov_tab$assignment_gene2) %in% c('LA', 'AL'), 'identity']
identities_of_LXs <- paralog_cov_tab[paste0(paralog_cov_tab$assignment_gene1, paralog_cov_tab$assignment_gene2) %in% c('LX', 'XL'), 'identity']
identities_of_AXs <- paralog_cov_tab[paste0(paralog_cov_tab$assignment_gene1, paralog_cov_tab$assignment_gene2) %in% c('AX', 'XA'), 'identity']

pal <- c('blue', 'yellow', 'green', 'orange', 'purple', 'red')

png('tables/paralog_divergence_hist.png')
hist(paralog_cov_tab$identity, breaks = 100, main = 'nt identity of genome-wide paralogs', xlab = 'Identity [ % ]')
hist(identities_of_LAs, col = 'blue', add = T, breaks = 100)
hist(identities_of_AAs, col = 'green', add = T, breaks = 100)
hist(identities_of_LLs, col = 'yellow', add = T, breaks = 50)
hist(identities_of_LXs, col = 'orange', add = T, breaks = 50)
hist(identities_of_AXs, col = 'purple', add = T, breaks = 50)
hist(identities_of_XXs, col = 'red', add = T, breaks = 50)
legend('topleft', bty = 'n', c('all', 'LA', 'LL', 'AA', 'LX', 'AX', 'XX'), col = c('black', pal), pch = 20)
dev.off()

len.LL<-length(identities_of_LLs)
len.AA<-length(identities_of_AAs)
len.XX<-length(identities_of_XXs)
len.LX<-length(identities_of_LXs)
len.AX<-length(identities_of_AXs)
len.LA<-length(identities_of_LAs)

len.paralogtypes<-c(len.AA, len.AX, len.LA, len.LL, len.LX, len.XX)
as.numeric(len.paralogtypes)# sum 2796 paralogs
paralog.types<-c('AA', 'AX', 'LA', 'LL', 'LX', 'XX')
paralog.nums<-data.frame(paralog.types, len.paralogtypes)
pal <- c('green', 'purple', 'blue', 'yellow', 'orange', 'red')
png('tables/paralog_count_all.png')
barplot(paralog.nums$len.paralogtypes~paralog.nums$paralog.types, col=pal)
dev.off()

#(1059+603+287)/(696+60+51+1059+603+287)= 70% of paralogs have something to do with L chromosome.
### Now some filtering of paralogs

#### filtering based on relative size of genes (only keeping genes that are withing 2X size of each other)
filter_paralog_cov_tab<-paralog_cov_tab[paralog_cov_tab$lendif1<=2,]
filter_paralog_cov_tab<-filter_paralog_cov_tab[filter_paralog_cov_tab$lendif2<=2,]
#dim 3399   18
filter_paralog_cov_tab<-filter_paralog_cov_tab[filter_paralog_cov_tab$percentaln1 >0.7 , ]
filter_paralog_cov_tab<-filter_paralog_cov_tab[filter_paralog_cov_tab$percentaln2 >0.7 , ]
#2707   18


ogtable3<-table(filter_paralog_cov_tab$og)
hist(ogtable3, breaks=30)
ogtable4<-as.data.frame(ogtable3)
filter_paralog_cov_tab<-merge(filter_paralog_cov_tab , ogtable4 ,by.x="og",by.y = "Var1")
head(filter_paralog_cov_tab)
colnames(filter_paralog_cov_tab)[21] <- 'og_freq2'

identities_of_LLs <- filter_paralog_cov_tab[paste0(filter_paralog_cov_tab$assignment_gene1, filter_paralog_cov_tab$assignment_gene2) == 'LL', 'identity']
identities_of_XXs <- filter_paralog_cov_tab[paste0(filter_paralog_cov_tab$assignment_gene1, filter_paralog_cov_tab$assignment_gene2) == 'XX', 'identity']
identities_of_AAs <- filter_paralog_cov_tab[paste0(filter_paralog_cov_tab$assignment_gene1, filter_paralog_cov_tab$assignment_gene2) == 'AA', 'identity']
identities_of_LAs <- filter_paralog_cov_tab[paste0(filter_paralog_cov_tab$assignment_gene1, filter_paralog_cov_tab$assignment_gene2) %in% c('LA', 'AL'), 'identity']
identities_of_LXs <- filter_paralog_cov_tab[paste0(filter_paralog_cov_tab$assignment_gene1, filter_paralog_cov_tab$assignment_gene2) %in% c('LX', 'XL'), 'identity']
identities_of_AXs <- filter_paralog_cov_tab[paste0(filter_paralog_cov_tab$assignment_gene1, filter_paralog_cov_tab$assignment_gene2) %in% c('AX', 'XA'), 'identity']


cbbPalette <- c("#E69F00", "#0072B2", "#009E73", "#F0E442","#56B4E9" , "#D55E00", "#CC79A7", "#000000")

filter_paralog_cov_tab$paralog <- paste0(filter_paralog_cov_tab$assignment_gene1, filter_paralog_cov_tab$assignment_gene2)
filter_paralog_cov_tab[filter_paralog_cov_tab$paralog=="LA","paralog"] <- "AL"
filter_paralog_cov_tab[filter_paralog_cov_tab$paralog=="XL","paralog"] <- "LX"
filter_paralog_cov_tab[filter_paralog_cov_tab$paralog=="AX","paralog"] <- "XA"
#write.table(filter_paralog_cov_tab, 'data/filtered_paralog_tab.tsv', quote = F, sep = "\t", row.names = F)


LAsubset <- filter_paralog_cov_tab[filter_paralog_cov_tab$paralog %in% c("AA", "AL", "LL"),]
cbbPalette <- c( "#009E73","#F0E442", "#E69F00","#CC79A7" ,"#56B4E9" , "#D55E00","#0072B2" , "#000000")
ggplot(LAsubset, aes(identity, fill = paralog)) + geom_histogram( binwidth=1)+ scale_fill_manual(values=cbbPalette) +
  theme_classic()+theme(axis.text=element_text(size=12))
  #call geom_histogram with position="dodge" to offset the bars and manual binwidth of 2
full.subset <- filter_paralog_cov_tab[filter_paralog_cov_tab$paralog %in% c("AA", "AL", "LL", "LX", "XX", "XA"),]
paralog.numbers<-table(full.subset$paralog)
paralog.numbers<-as.data.frame(paralog.numbers)


## plots for paper (fig3)
cbbPalette <- c( "#009E73","#F0E442", "#E69F00","#CC79A7"  ,"#0072B2","#56B4E9", "#D55E00","#0072B2" , "#000000")
ggplot(paralog.numbers, aes(y=Freq, x=Var1, fill=Var1)) + scale_fill_manual(values=cbbPalette) +
  geom_bar(position="dodge", stat="identity")+theme_classic()+theme(axis.text=element_text(size=14))
ggplot(full.subset, aes(identity, fill = paralog)) + geom_histogram( binwidth=1) + scale_fill_manual(values=cbbPalette) +
  theme_classic()+theme(axis.text=element_text(size=12))
#+ theme(legend.text=element_text(size=12))
ggplot(full.subset, aes(identity, fill = paralog)) + geom_histogram( binwidth=1)+ xlim(70,100)+ scale_fill_manual(values=cbbPalette) +
  theme_classic()+theme(axis.text=element_text(size=12))


table(LAsubset$paralog)
 


pal <- c('blue', 'yellow', 'green', 'orange', 'purple', 'red')

png('tables/paralog_divergence_hist_filtered.png')
hist(paralog_cov_tab$identity, breaks = 100, main = 'nt identity of genome-wide paralogs', xlab = 'Identity [ % ]')
hist(identities_of_LAs, col = 'blue', add = T, breaks = 100)
hist(identities_of_AAs, col = 'green', add = T, breaks = 50)
hist(identities_of_LLs, col = 'yellow', add = T, breaks = 50)
hist(identities_of_LXs, col = 'orange', add = T, breaks = 50)
hist(identities_of_AXs, col = 'purple', add = T, breaks = 50)
hist(identities_of_XXs, col = 'red', add = T, breaks = 50)
legend('topleft', bty = 'n', c('all', 'LA', 'LL', 'AA', 'LX', 'AX', 'XX'), col = c('black', pal), pch = 20)
dev.off()

len.LL<-length(identities_of_LLs)
len.AA<-length(identities_of_AAs)
len.XX<-length(identities_of_XXs)
len.LX<-length(identities_of_LXs)
len.AX<-length(identities_of_AXs)
len.LA<-length(identities_of_LAs)

len.paralogtypes<-c(len.AA, len.AX, len.LA, len.LL, len.LX, len.XX)
as.numeric(len.paralogtypes)# sum 2796 paralogs
paralog.types<-c('AA', 'AX', 'LA', 'LL', 'LX', 'XX')
paralog.nums<-data.frame(paralog.types, len.paralogtypes)
pal <- c('green', 'purple', 'blue', 'yellow', 'orange', 'red')
png('tables/paralog_count_all_filtered.png')
barplot(paralog.nums$len.paralogtypes~paralog.nums$paralog.types, col=pal)
dev.off()


### Looking at all L-L paralogs (in filtered list)
#1. to look at how many L-L paralogs there are with no other paralogs
#2. to look at whether there is evidence for two distinct L chromosomes (maybe could also look at how many non L_L paralogs have a og_freq>1)


#Now I want to see how many it takes out if I get rid of the paralogs that hit to dif genes on the same contig

#dim filtered table= 2707   19
filter_paralog_cov_tab_difscf<-filter_paralog_cov_tab[filter_paralog_cov_tab$scf_gene1!=filter_paralog_cov_tab$scf_gene2, ]
dim(filter_paralog_cov_tab_difscf)
#dim2469   19
identities_of_LLs <- filter_paralog_cov_tab_difscf[paste0(filter_paralog_cov_tab_difscf$assignment_gene1, filter_paralog_cov_tab_difscf$assignment_gene2) == 'LL', 'identity']
identities_of_XXs <- filter_paralog_cov_tab_difscf[paste0(filter_paralog_cov_tab_difscf$assignment_gene1, filter_paralog_cov_tab_difscf$assignment_gene2) == 'XX', 'identity']
identities_of_AAs <- filter_paralog_cov_tab_difscf[paste0(filter_paralog_cov_tab_difscf$assignment_gene1, filter_paralog_cov_tab_difscf$assignment_gene2) == 'AA', 'identity']
identities_of_LAs <- filter_paralog_cov_tab_difscf[paste0(filter_paralog_cov_tab_difscf$assignment_gene1, filter_paralog_cov_tab_difscf$assignment_gene2) %in% c('LA', 'AL'), 'identity']
identities_of_LXs <- filter_paralog_cov_tab_difscf[paste0(filter_paralog_cov_tab_difscf$assignment_gene1, filter_paralog_cov_tab_difscf$assignment_gene2) %in% c('LX', 'XL'), 'identity']
identities_of_AXs <- filter_paralog_cov_tab_difscf[paste0(filter_paralog_cov_tab_difscf$assignment_gene1, filter_paralog_cov_tab_difscf$assignment_gene2) %in% c('AX', 'XA'), 'identity']

pal <- c('blue', 'yellow', 'green', 'orange', 'purple', 'red')

png('tables/paralog_divergence_hist_filtered_difscf.png')
hist(paralog_cov_tab$identity, breaks = 100, main = 'nt identity of genome-wide paralogs', xlab = 'Identity [ % ]')
hist(identities_of_LAs, col = 'blue', add = T, breaks = 100)
hist(identities_of_AAs, col = 'green', add = T, breaks = 100)
hist(identities_of_LLs, col = 'yellow', add = T, breaks = 50)
hist(identities_of_LXs, col = 'orange', add = T, breaks = 50)
hist(identities_of_AXs, col = 'purple', add = T, breaks = 50)
hist(identities_of_XXs, col = 'red', add = T, breaks = 50)
legend('topleft', bty = 'n', c('all', 'LA', 'LL', 'AA', 'LX', 'AX', 'XX'), col = c('black', pal), pch = 20)
dev.off()

len.LL<-length(identities_of_LLs)
len.AA<-length(identities_of_AAs)
len.XX<-length(identities_of_XXs)
len.LX<-length(identities_of_LXs)
len.AX<-length(identities_of_AXs)
len.LA<-length(identities_of_LAs)

len.paralogtypes<-c(len.AA, len.AX, len.LA, len.LL, len.LX, len.XX)
as.numeric(len.paralogtypes)# sum 2796 paralogs
paralog.types<-c('AA', 'AX', 'LA', 'LL', 'LX', 'XX')
paralog.nums<-data.frame(paralog.types, len.paralogtypes)
pal <- c('green', 'purple', 'blue', 'yellow', 'orange', 'red')
png('tables/paralog_count_all_filtered_difscf.png')
barplot(paralog.nums$len.paralogtypes~paralog.nums$paralog.types, col=pal)
dev.off()
#filter_paralog_cov_tab_difscf[filter_paralog_cov_tab_difscf$scf_gene1==filter_paralog_cov_tab_difscf$scf_gene2, ]




LL_paralogs<-filter_paralog_cov_tab_difscf[filter_paralog_cov_tab_difscf$assignment_gene1=="L", ]
LL_paralogs<-LL_paralogs[LL_paralogs$assignment_gene2=="L", ]
LL_paralogs <- LL_paralogs[!is.na(LL_paralogs$identity),]
#not sure why but this filtering is for some reason adding rows of NA's into the data frame
#dim(LL_paralogs) # 382  19

lowcov_L <- apply(LL_paralogs[,c(17,18)], 1, min)
highcov_L <- apply(LL_paralogs[,c(17,18)], 1, max)

#and plot the two histograms
cbbPalette <- c("#E69F00", "#0072B2", "#009E73", "#F0E442","#56B4E9" , "#D55E00", "#CC79A7", "#000000")

hist(highcov_L, col = c("#999999"), breaks = 1000, xlim=c(0,80), ylim = c(0,50))
hist(lowcov_L, col = rgb(0.90,0.62,0, alpha = 0.7), breaks=100, xlim=c(0,100),add = T)
hist(highcov_L- lowcov_L, breaks = 2000, xlim = c(0,30))
###note: there are a few with a cov greater than 200, look at those

table(LL_paralogs$og_freq1)
table(LL_paralogs$og_freq2)
#looking at genes with really high cov levels
#LL_paralogs[LL_paralogs$meancov_gene1>100, ]

##below here I'm filtering only L genes with a reasonable cov based on assembly <50, and looking at the hist of these pairs
LL_paralogs_cov50<-LL_paralogs[LL_paralogs$meancov_gene1<50, ]
LL_paralogs_cov50<-LL_paralogs_cov50[LL_paralogs_cov50$meancov_gene2<50, ]
#> dim(LL_paralogs_cov50)[1] 328  20

lowcov_L <- apply(LL_paralogs_cov50[,c(17,18)], 1, min)
highcov_L <- apply(LL_paralogs_cov50[,c(17,18)], 1, max)

cbbPalette <- c("#E69F00", "#0072B2")
hist(highcov_L, col = c("#D55E00"), breaks = 30, xlim=c(0,60), ylim = c(0,50), xlab="Gene Coverage", cex.lab=1.2)
hist(lowcov_L, col = rgb(0.90,0.62,0, alpha = 0.65), breaks=30, xlim=c(0,100),add = T)

#hist(highcov_L, col = rgb(0.80,0.40,0, alpha = 0.8), breaks = 30, xlim=c(0,60), ylim = c(0,50))
#hist(lowcov_L, col = rgb(0.95,0.9,0.25, alpha = 0.7), breaks=30, xlim=c(0,100),add = T)


#> mean(highcov_L)[1] 30.34548, > mean(lowcov_L)[1] 25.38722
#> sd(highcov_L)[1] 3.752866,   > sd(lowcov_L)[1] 3.576065
#problem is that there needs to be a reason to set the cov values to a certain amount max







### ok, so I've decided for now that I'm going to filter based on the relative length of genes involved in the blast hit (2X) and the proportion of each gene that the alignment covers (70%)
#### and then I'll set the og group freq's afterwards (since the point of filtering is only to get reliable alignments)

#### now I'm going to try to figure out how many genes are involved in each OG group and how that compares with the og freq (i.e. will give an indication of whether there are some genes in the group that don't blast to each other)

#filter_paralog_cov_tab
#filtering above paralogs so we only have those assigned as L, X, or A (not c or NA ones)
#final_filtered_paralogs <- filter_paralog_cov_tab[(filter_paralog_cov_tab$assignment_gene1, filter_paralog_cov_tab$assignment_gene2) %in% c('LA', 'AL'), 'identity']
final_filtered_paralogs <- filter_paralog_cov_tab[ which( filter_paralog_cov_tab$assignment_gene1=="A" | filter_paralog_cov_tab$assignment_gene1=="L" | filter_paralog_cov_tab$assignment_gene1=="X") , ]
final_filtered_paralogs <- final_filtered_paralogs[ which( final_filtered_paralogs$assignment_gene2=="A" | final_filtered_paralogs$assignment_gene2=="L" | final_filtered_paralogs$assignment_gene2=="X") , ]
#dim(final_filtered_paralogs)2210   21
#dim(filter_paralog_cov_tab)2707   21
final_filtered_paralogs_difscf<-final_filtered_paralogs[final_filtered_paralogs$scf_gene1!=final_filtered_paralogs$scf_gene2, ]
#> dim(final_filtered_paralogs_difscf)[1] 1983   20

gene.og<-final_filtered_paralogs_difscf[, c(1,2)]
gene.og2<-final_filtered_paralogs_difscf[, c(1,3)]
colnames(gene.og)[2] <- 'gene'
colnames(gene.og2)[2] <- 'gene'
#getting the number of genes in each og reciprocal blast search (some can be there twice)
gene.og<-rbind(gene.og,gene.og2)
og.gene.num<-table(gene.og$og)
og.gene.num<-as.data.frame(og.gene.num)
colnames(og.gene.num)[1] <- 'og'
colnames(og.gene.num)[2] <- 'total.genes'
#getting the number of genes in each og group (each gene can only be there once)
gene.og1<-unique(gene.og)
og.gene.num1<-table(gene.og1$og)
og.gene.num1<-as.data.frame(og.gene.num1)
colnames(og.gene.num1)[1] <- 'og'
colnames(og.gene.num1)[2] <- 'unique.genes'

og.group.summary<-merge(og.gene.num , og.gene.num1 ,by.x="og",by.y = "og")

##filtering so I get only total genes=6 and unique genes=3 
#(I think that should give me just perfect pairs with 6 genes)
paralog.3gene<-og.group.summary[og.group.summary$total.genes==6,]
paralog.3gene1<-paralog.3gene[paralog.3gene$unique.genes==3,]
##ok, so this is the same answer I got the first time but this time at least the code works again
##code works until here, don't run below
#If I want more info on paralog with 3 genes I need to pull them out of the full list
para.3pairs<-merge(paralog.3gene1, final_filtered_paralogs, by.x="og",by.y = "og")
##ok, there is a problem with this list, it includes all NA and c assignments, need to take them out first)
##also need to take out ones where the node for the query and subject are the same
#> dim(para.3pairs) 234  23
#> dim(para.3pairs) [1] 183  22 (after filtering paralogs on same contig)
para.3pairs<-para.3pairs[, c(1:13,19,20)]

#histogram of types of paralogs

