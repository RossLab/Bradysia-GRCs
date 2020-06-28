# Now we have some all-by-all reciprocal blast results. We want to add some extra info to the blast results table
### there are several blast outputs that we generated
### 1. reciprocal blast hits table with each recirocal pair of genes, similarity, and length of alignment, and gene lengths
### 2. paralog groups with all reciprocal blast hits that contains all the genes that are paralogs to each other
#### we did this for aa and nt sequences and have tables at different similarity thresholds for the aa blasts

## BELOW I'm looking at just the nt blast results
#### note: the blast parameters were set so only hits with an e-value < 1e-10 are kept

### First making a table with gene ident, length, assignment, etc...(need to get gene info from braker list)

#setwd("/Users/christina/projects/Sciara-L-chromosome/")
scf_asn <- read.delim("data/scaffold_assignment_tab_full.tsv", header=T, stringsAsFactors = F)
trans_asn <- read.delim("data/gene.scaffold.map.tsv", header=F, stringsAsFactors = F, col.names = c('gene', 'scf'))
row.names(scf_asn) <- scf_asn$scf
row.names(trans_asn) <- trans_asn$gene
genes.to.scf<-merge(scf_asn,trans_asn)

genes.info<-genes.to.scf[,c(16,15,1)]

### now loading some of the blast results and merging data.frames into one table

#### nt_blast results
blast_nt <- read.delim("data/nt_orthology_OG_pairs.tsv", header=T, stringsAsFactors = F, col.names = c("og", "gene1" , "gene2", "identity", "aln_length", "len_gene1", "len_gene2"))
paralog_table_full<- merge(blast_nt ,genes.info,by.x="gene1",by.y = "gene")
#paralog_table_full<-paralog_table_full[,c(1:9)]
colnames(paralog_table_full)[8] <- 'assignment_gene1'
colnames(paralog_table_full)[9] <- 'scf_gene1'

paralog_table_full<- merge(paralog_table_full ,genes.info,by.x="gene2",by.y = "gene")
#paralog_table_full<-paralog_table_full[,c(1:11)]
colnames(paralog_table_full)[10] <- 'assignment_gene2'
colnames(paralog_table_full)[11] <- 'scf_gene2'

#paralog_table_full<-paralog_table_full[, c(3,2,1,4,5,6,8,9,7,10,11)]

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
colnames(paralog_cov_tab)[17] <- 'meancov_gene1'
paralog_cov_tab <- merge(paralog_cov_tab,cov_genes,by.x="gene2",by.y = "braker_gene")
colnames(paralog_cov_tab)[18] <- 'meancov_gene2'
colnames(paralog_cov_tab)[16] <- 'og_freq1'
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

#### I think I might also want to split up the lists based on frequeny (which is the number of paralogs), I'm not sure it I should do 3 and less or only 1 and more than one)

table(paralog_cov_tab$Freq)
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
colnames(filter_paralog_cov_tab)[19] <- 'og_freq2'

identities_of_LLs <- filter_paralog_cov_tab[paste0(filter_paralog_cov_tab$assignment_gene1, filter_paralog_cov_tab$assignment_gene2) == 'LL', 'identity']
identities_of_XXs <- filter_paralog_cov_tab[paste0(filter_paralog_cov_tab$assignment_gene1, filter_paralog_cov_tab$assignment_gene2) == 'XX', 'identity']
identities_of_AAs <- filter_paralog_cov_tab[paste0(filter_paralog_cov_tab$assignment_gene1, filter_paralog_cov_tab$assignment_gene2) == 'AA', 'identity']
identities_of_LAs <- filter_paralog_cov_tab[paste0(filter_paralog_cov_tab$assignment_gene1, filter_paralog_cov_tab$assignment_gene2) %in% c('LA', 'AL'), 'identity']
identities_of_LXs <- filter_paralog_cov_tab[paste0(filter_paralog_cov_tab$assignment_gene1, filter_paralog_cov_tab$assignment_gene2) %in% c('LX', 'XL'), 'identity']
identities_of_AXs <- filter_paralog_cov_tab[paste0(filter_paralog_cov_tab$assignment_gene1, filter_paralog_cov_tab$assignment_gene2) %in% c('AX', 'XA'), 'identity']

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

### Looking at all L-L paralogs (in filtered list)
#1. to look at how many L-L paralogs there are with no other paralogs
#2. to look at whether there is evidence for two distinct L chromosomes (maybe could also look at how many non L_L paralogs have a og_freq>1)

LL_paralogs<-filter_paralog_cov_tab_difscf[filter_paralog_cov_tab_difscf$assignment_gene1=="L", ]
LL_paralogs<-LL_paralogs[LL_paralogs$assignment_gene2=="L", ]
LL_paralogs <- LL_paralogs[!is.na(LL_paralogs$identity),]
#not sure why but this filtering is for some reason adding rows of NA's into the data frame
#dim(LL_paralogs) # 382  19

lowcov_L <- apply(LL_paralogs[,c(17,18)], 1, min)
highcov_L <- apply(LL_paralogs[,c(17,18)], 1, max)

#and plot the two histograms

hist(highcov_L, col = 'purple', breaks = 1000, xlim=c(0,100), ylim = c(0,50))
hist(lowcov_L, col = rgb(0.4,0.6,0.1, alpha = 0.8), breaks=100, xlim=c(0,100),add = T)
hist(highcov_L- lowcov_L, breaks = 2000, xlim = c(0,30))
###note: there are a few with a cov greater than 200, look at those

table(LL_paralogs$og_freq1)
table(LL_paralogs$og_freq2)
#looking at genes with really high cov levels
#LL_paralogs[LL_paralogs$meancov_gene1>100, ]

### ok, so I've decided for now that I'm going to filter based on the relative length of genes involved in the blast hit (2X) and the proportion of each gene that the alignment covers (70%)
#### and then I'll set the og group freq's afterwards (since the point of filtering is only to get reliable alignments)

#### now I'm going to try to figure out how many genes are involved in each OG group and how that compares with the og freq (i.e. will give an indication of whether there are some genes in the group that don't blast to each other)

filter_paralog_cov_tab

gene.og<-filter_paralog_cov_tab[, c(1,2)]
gene.og2<-filter_paralog_cov_tab[, c(1,3)]

colnames(gene.og)[2] <- 'gene'
colnames(gene.og2)[2] <- 'gene'
gene.og<-rbind(gene.og,gene.og2)
gene.og<-unique(gene.og)
og.gene.num<-table(gene.og$og)
og.gene.num<-as.data.frame(og.gene.num)
filter_paralog_cov_tab<-merge(filter_paralog_cov_tab , og.gene.num ,by.x="og",by.y = "Var1")
colnames(filter_paralog_cov_tab)[20]<-'genes.paralog.group'

og.group.filt<-filter_paralog_cov_tab[, c(1,19,20)]
colnames(og.group.filt)[3] <- 'genes'
og.group.filt$expected_links = sapply(og.group.filt$genes, function(x) { ncol(combn(x, 2))} )
og.group.filt1<-og.group.filt[og.group.filt$og_freq2 == og.group.filt$expected_links, ]
table(og.group.filt1$genes) 
perfect.og<-unique(og.group.filt1)
perfect.og.tab<-merge(filter_paralog_cov_tab , perfect.og ,by.x="og",by.y = "og")

########### I want to see whether the numbers are a lot different if I don't filter the hits (i.e. if I'm losing perfect paralog groups by filtering)

gene.og_unf<-paralog_cov_tab[, c(1,3)]
gene.og2_unf<-paralog_cov_tab[, c(2,3)]

colnames(gene.og_unf)[1] <- 'gene'
colnames(gene.og2_unf)[1] <- 'gene'
gene.og_unf<-rbind(gene.og_unf,gene.og2_unf)
gene.og_unf<-unique(gene.og_unf)
og.gene.num_unf<-table(gene.og_unf$og)
og.gene.num_unf<-as.data.frame(og.gene.num_unf)
paralog_cov_tab_unf<-merge(paralog_cov_tab , og.gene.num_unf ,by.x="og",by.y = "Var1")


og.group.filt_unf<-paralog_cov_tab_unf[, c(1,16,19)]
colnames(og.group.filt_unf)[3] <- 'genes'
og.group.filt_unf$expected_links = sapply(og.group.filt_unf$genes, function(x) { ncol(combn(x, 2))} )
og.group.filt_unf<-og.group.filt_unf[og.group.filt_unf$og_freq1 == og.group.filt_unf$expected_links, ]
table(og.group.filt_unf$genes)  

#doing some data exploration
paralog.3gene<-perfect.og.tab[perfect.og.tab$genes==3,]

#identities_of_AAAs <- filter_paralog_cov_tab_difscf[paste0(filter_paralog_cov_tab_difscf$assignment_gene1, filter_paralog_cov_tab_difscf$assignment_gene2) %in% c('AX', 'XA'), 'identity']


   
