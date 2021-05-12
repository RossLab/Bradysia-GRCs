#setwd("/Users/christina/projects/Sciara-L-chromosome/")
library('ggplot2')
ax_genes <- read.table('tables/A-busco-phylogenies-summary.tsv', header = T)
busco_asn<- read.table('tables/BUSCO_assigned.Illumina.tsv', header = T, sep = '\t', stringsAsFactors = F, quote = "\"")

# Histogram of scp AX BUSCO placement in phylogenies
ax_genes<-ax_genes[c(-126),]
phy_asn<-as.data.frame(table(ax_genes$Core_genome))
colnames(phy_asn)[1]<-"clade"
colnames(phy_asn)[2]<-"num.genes"
phy_asn$clade<-factor(phy_asn$clade, levels=c("sciaridae", "cecidomyiidae", "other"))
#cecidomyiidae=1         other=22     sciaridae=585

pu <- rgb(206,146,218, maxColorValue = 255)#ceci
te <- rgb(76,206,175, maxColorValue = 255)#sciaridae
#bl <- rgb(44,127,184, maxColorValue = 255)
pal <- c(te, pu, 'grey')

pdf('figures/phylo_trees_topologies_AX.pdf', width = 2, height = 4)
barplot(phy_asn$num.genes~phy_asn$clade, col = pal, cex.axis = 0.8, ylim=c(0,600))
dev.off()

# branch distribution differences between A and X scp BUSCO genes
busco_asn<-busco_asn[,c(1,4)]
branch_asn <- merge(ax_genes,busco_asn,by.x="BUSCO_id",by.y = "busco_id")
branch_asn <- branch_asn[branch_asn$Core_genome=="sciaridae",]
boxplot(branch_asn$Core_len~branch_asn$assignment)
#206,146,218
#76,206,175
greenA <- rgb(0, 158, 115, maxColorValue = 255)
blueX <- rgb(0, 114, 178, maxColorValue = 255)

light_greenA <- rgb(t(col2rgb(greenA)), maxColorValue=255, alpha = 150)
light_blueX <- rgb(t(col2rgb(blueX)), maxColorValue=255, alpha = 150)


A_branch_lengths <- c(branch_asn[branch_asn$assignment == "A", "Core_len"])
X_branch_lengths <- c(branch_asn[branch_asn$assignment == "X", "Core_len"])

df <- data.frame(branch_length = c(A_branch_lengths, X_branch_lengths),
                 branch_type = c(rep('A', length(A_branch_lengths)), rep('X', length(X_branch_lengths))))
mu <- data.frame(branch_type = c('A', 'X'), grp.mean = c(mean(A_branch_lengths), mean(X_branch_lengths)))
p <- ggplot(df, aes(x=branch_length, color=branch_type, fill = branch_type)) +
  geom_density() +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=branch_type),
             linetype="dashed")


pdf('figures/AXscp_busco_branch_length_distribution.pdf', width = 7, height = 4)

p + scale_color_manual(values = c(light_greenA, light_blueX)) + scale_fill_manual(values=c(light_greenA, light_blueX)) +theme_classic()

dev.off()

wilcox.test(A_branch_lengths,X_branch_lengths)
