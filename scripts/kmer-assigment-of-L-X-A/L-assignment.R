library(ggplot2)
# mapped kmers
#setwd("/Users/christina/projects/Sciara-L-chromosome/")
sr_kmers <- read.table("data/mapping/table_of_mapped_kmers_spades.tsv", header = T, sep = '\t', stringsAsFactors = F)
colnames(sr_kmers)[1] <- "ID"

coverage_tab <- read.table('data/table.covdiff.germ.soma.txt', stringsAsFactors = F, header = T)
coverage_tab <- coverage_tab[coverage_tab$ID != "*",]
colnames(coverage_tab)[2] <- 'len'

coverage_tab$scf <- unlist(lapply(strsplit(coverage_tab$ID, "_"), function(x){ paste(x[1:2], collapse = "_") } ))
coverage_tab <- coverage_tab[,c(6, 1:5)]

scf_to_add <- coverage_tab$ID[!coverage_tab$ID %in% sr_kmers$ID]
number_to_add <- length(scf_to_add)
sr_kmers <- rbind(sr_kmers, data.frame(ID = scf_to_add,
                                       L = rep(0, number_to_add),
                                       X = rep(0, number_to_add),
                                       A = rep(0, number_to_add)))

overlap_tab <- merge(coverage_tab, sr_kmers)
# row.names(overlap_tab) <- overlap_tab$ID

overlap_tab$kmer_naive <- apply(overlap_tab[,c('L','X','A')], 1, function(x){ c('L', 'X', 'A')[which.max(x)] })
overlap_tab$score <- apply(overlap_tab[,c('L','X','A')], 1, function(x){ max(x) / sum(x) })
overlap_tab$L_score <- overlap_tab$L / overlap_tab$len
overlap_tab$A_score <- overlap_tab$A / overlap_tab$len
overlap_tab$X_score <- overlap_tab$X / overlap_tab$len

overlap_tab <- overlap_tab[order(overlap_tab$len, decreasing = T),]

# overlap_tab$logdif2
L_assignment = overlap_tab$L_score > 0.8 & overlap_tab$logdif2 > 0.5
Lc_assignment = (overlap_tab$L_score > 0.8 | overlap_tab$logdif2 > 0.5) & !L_assignment
A_assignment = overlap_tab$A_score > 0.4 & overlap_tab$logdif2 < -0.1 & overlap_tab$logdif2 > -1
Ac_assignment = (overlap_tab$A_score > 0.4 | (overlap_tab$logdif2 < -0.1 & overlap_tab$logdif2 > -1)) & !A_assignment
X_assignment = overlap_tab$X_score > 0.4 & overlap_tab$logdif2 < 0.5 & overlap_tab$logdif2 > -0.1
Xc_assignment = (overlap_tab$X_score > 0.4 | (overlap_tab$logdif2 < 0.5 & overlap_tab$logdif2 > -0.1)) & !X_assignment

categories <- c('L', 'Lc', 'A', 'Ac', 'X', 'Xc', 'NA')
overlap_tab$assignments <- 'NA'
overlap_tab$assignments[L_assignment] <- 'L'
overlap_tab$assignments[Lc_assignment] <- 'Lc'
overlap_tab$assignments[X_assignment] <- 'X'
overlap_tab$assignments[Xc_assignment] <- 'Xc'
overlap_tab$assignments[A_assignment] <- 'A'
overlap_tab$assignments[Ac_assignment] <- 'Ac'
overlap_tab$assignments[(Ac_assignment + Xc_assignment + Lc_assignment) > 1] <- 'NA'

#write.table(overlap_tab, 'data/scaffold_assignment_tab_full.tsv', quote = F, sep = "\t", row.names = F)


# L_kmer_candidates <- overlap_tab[L_assignment, c('ID', 'len', 'L_score', 'logdif2')]
# coverage_candidates <- overlap_tab[cov_assignment, c('ID', 'len', 'L_score', 'logdif2')]
#
# shared_candidates <- merge(coverage_candidates, L_kmer_candidates)
# kmer_specific_candidates <- L_kmer_candidates[!L_kmer_candidates$ID %in% coverage_candidates$ID,]
# coverage_specific_candidates <- coverage_candidates[!coverage_candidates$ID %in% L_kmer_candidates$ID,]

"table of assignments [Mbp]"

sizes <- round(sapply(categories, function(x){ sum(overlap_tab[overlap_tab$assignments == x, 'len']) / 1e6 } ), 1)

cat(paste0('|    ', paste(categories, collapse = '    |    '), '    |\n'))
cat(paste0('| ---', paste(sapply(categories, function(x){ paste(rep("-", nchar(x)), collapse = "") }), collapse = '--- | ---'), '--- |\n'))
cat(paste0('|', paste(sapply(sizes, function(x){ paste(c(" ", x, rep(" ", 8 - nchar(x))), collapse = "") }), collapse = '| '), ' |\n'))


#I want to make a small plot of how many of each thing we have (i.e. total length)

pal <- c('yellow','yellow', 'green','green', 'red', 'red', 'grey')
barplot(sizes, ylab = "size (Mb)", xlab='chromosome assignment',col=pal)

# | First Header  | Second Header |
# | ------------- | ------------- |
# | Content Cell  | Content Cell  |
# | Content Cell  | Content Cell  |

# "all putative L [Mbp]"
# sum(overlap_tab[overlap_tab$assignments %in% c('L', 'Lk', 'Lc', 'Lp'), 'len']) / 1e6
# "both [Mbp]"
# sum(shared_candidates$len) / 1e6
# "kmers only [Mbp]"
# sum(kmer_specific_candidates$len) / 1e6
# "coverage only [Mbp]"
# sum(coverage_specific_candidates$len) / 1e6
# "all putative A [Mbp]"

# shared_candidates <- shared_candidates[order(shared_candidates$len, decreasing = T),]
# kmer_specific_candidates <- kmer_specific_candidates[order(kmer_specific_candidates$len, decreasing = T),]
# coverage_specific_candidates <- coverage_specific_candidates[order(coverage_specific_candidates$len, decreasing = T),]
#
# write.table(shared_candidates[,'ID'], file = "data/L-candidates-both.tsv",
#             quote = F, sep = "\t", row.names = F, col.names = F)
# write.table(kmer_specific_candidates[,'ID'], file = "data/L-candidates-kmers-only.tsv",
#             quote = F, sep = "\t", row.names = F, col.names = F)
# write.table(coverage_specific_candidates[,'ID'], file = "data/L-candidates-cov-only.tsv",
#             quote = F, sep = "\t", row.names = F, col.names = F)

### I want to make a summary table comparing chromosome sizes between Rasch and our assembly and identification
sizes<- as.data.frame(sizes)
sizes$identity<-c("L", "Lc", "A", "Ac", "X", "Xc", "NA")
sizes2<-sizes[c(1,3,5),]
6.8+0.9+2.4+10.1
#rbind<-(sizes2 , na=c(20.2, "unassigned"))

chromosome.type<-c("A","X", "L","Unassigned", "Total","A","X", "L","Unassigned", "Total")
size<-c(162.4, 52.9, 154.1, 20.2, 389.6 ,225, 49, 88, 0, 362)
source<-c("genome","genome","genome","genome","genome","rasch","rasch","rasch","rasch","rasch")

size.table<-cbind(chromosome.type,size.assembly,size.rasch)
pal <- c('yellow','yellow', 'green','green', 'red', 'red', 'grey')
# The palette with black:
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
barplot(size.assembly, ylab = "size (Mb)", xlab='chromosome assignment',col=cbbPalette)

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2","#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
chromosome<-c("Autosome","X", "GRC","Unassigned", "Total","Austosome","X", "GRC","Unassigned", "Total")
tab.size<-as.data.frame(chromosome)
tab.size$size<-c(162.4, 52.9, 154.1, 20.2, 389.6 ,225, 49, 88, 0, 362)
tab.size$source<-c("genome","genome","genome","genome","genome","rasch","rasch","rasch","rasch","rasch")

ggplot(tab.size, aes(fill=chromosome, y=size, x=source)) + scale_fill_manual(values=cbbPalette) +
  geom_bar(position="stack", stat="identity")+theme_classic()

ggplot(tab.size, aes(fill=chromosome, y=size, x=source)) + scale_fill_manual(values=cbbPalette) +
  geom_bar(position="dodge", stat="identity")+theme_classic()


#without the total bar
cbbPalette <- c("#009E73", "#56B4E9","#E69F00","#999999", "#F0E442","#0072B2" , "#D55E00", "#CC79A7")
chromosome<-c("Autosomes","X", "GRC","Unassigned","Autosomes","X", "GRC","Unassigned")
tab.size<-as.data.frame(chromosome)
tab.size$chromosome <- factor(tab.size$chromosome,levels = c("Autosomes", "X", "GRC", "Unassigned"))  
tab.size$size<-c(162.4, 52.9, 154.1, 20.2 ,225, 49, 88, 0)
tab.size$source<-c("genome","genome","genome","genome","rasch","rasch","rasch","rasch")
tab.size$source <- factor(tab.size$source,levels = c("rasch", "genome"))  

ggplot(tab.size, aes(fill=chromosome, y=size, x=source)) + scale_fill_manual(values=cbbPalette) +
  geom_bar(position="dodge", stat="identity")+theme_classic()+theme(axis.text=element_text(size=12))+ theme(legend.text=element_text(size=12))


cbbPalette <- c("#009E73", "#56B4E9","#E69F00","#999999", "#F0E442","#0072B2" , "#D55E00", "#CC79A7")
chromosome<-c("Autosomes","X", "GRC","Unassigned")
tab.size<-as.data.frame(chromosome)
tab.size$chromosome <- factor(tab.size$chromosome,levels = c("Autosomes", "X", "GRC", "Unassigned"))  
tab.size$size<-c(162.4, 52.9, 154.1, 20.2)
#tab.size$source<-c("genome","genome","genome","genome","rasch","rasch","rasch","rasch")
#tab.size$source <- factor(tab.size$source,levels = c("rasch", "genome"))  

ggplot(tab.size, aes(fill=chromosome, y=size, x=chromosome)) + scale_fill_manual(values=cbbPalette) +
  geom_bar(position="dodge", stat="identity")+theme_classic()+theme(axis.text=element_text(size=12))+ theme(legend.text=element_text(size=12))

