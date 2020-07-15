### Exploring colinearity data
library(ggplot2)

#setwd("/Users/christina/projects/Sciara-L-chromosome/")
colinear_tab <- read.table('data/collinear_genes_full_table.tsv', stringsAsFactors = F, header = T)

gene_info <- read.table('data/gene_cov_table.tsv', stringsAsFactors = F, header = T)
gene_info<- gene_info[, c(1,10,11)]

colinear_tab$paralog <- paste0(colinear_tab$chr, colinear_tab$paralog_asn)
colinear_tab[colinear_tab$paralog=="LA","paralog"] <- "AL"
colinear_tab[colinear_tab$paralog=="LX","paralog"] <- "XL"
colinear_tab<-merge(colinear_tab , gene_info ,by.x="gene",by.y = "braker_gene")
colinear_tab<-colinear_tab[,c(-7)]
#table(colinear_tab$paralog)
#para_typeAL=   LL   XL 
#num_para*2= 1440  484  288 
#num_blocks= 79  33  15


a_colinear_tab<-colinear_tab[grepl("*a", colinear_tab$block),]
b_colinear_tab<-colinear_tab[grepl("*b", colinear_tab$block),]

a_colinear_tab$block_id<-gsub('.{1}$', '', a_colinear_tab$block)
b_colinear_tab$block_id<-gsub('.{1}$', '', b_colinear_tab$block)
a_colinear_tab$side<-"a"
b_colinear_tab$side<-"b"
full_colinear_tab<-rbind(a_colinear_tab, b_colinear_tab)


#hist(LLa_colinear_tab$mean_cov, breaks=30)
#hist(LLb_colinear_tab$mean_cov, breaks=40, add=T)

LL_colinear_tab<-full_colinear_tab[full_colinear_tab$paralog=="LL",]
XL_colinear_tab<-full_colinear_tab[full_colinear_tab$paralog=="XL",]
AL_colinear_tab<-full_colinear_tab[full_colinear_tab$paralog=="AL",]

#scf_uniq<-LL_colinear_tab[!duplicated(LL_colinear_tab$orig_scf), ]
#scf_uniq<-scf_uniq[order(scf_uniq$block_id),]

#table(LL_colinear_tab$block_id)
#100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127  95  96  97  98  99 
# 12  12  22  16  12  18  16  12  14  14  16  12  12  12  12  22  20  12  18  12  12  16  12  18  16  22  14  12  14  12  16  12  12 
pal<-c("#999999", "#E69F00")
palette(pal)
for( block in unique(LL_colinear_tab$block_id)){
  pdf(paste('figures/LL_block_cov/', block ,'.pdf'))
  block_data <-LL_colinear_tab[LL_colinear_tab$block_id==block,]
  block_data$side <- as.factor(block_data$side)
  plot(mean_cov ~ order_in_block, data = block_data, pch = 20, col = side)
  dev.off()
}



#block96<-LL_colinear_tab[LL_colinear_tab$block_id=="96",]
#qplot(order_in_block, mean_cov, data=block96, colour=side)


#scf_uniq<-LL_colinear_tab[!duplicated(LL_colinear_tab$orig_scf), ]
#scf_uniq<-scf_uniq[order(scf_uniq$block_id),]
#scf_uniq<-scf_uniq[,c(3,4,9,12,13)]
#write.table(scf_uniq, 'data/colinear_scf.tsv', quote = F, sep = "\t", row.names = F)

#write.table(full_colinear_tab, 'data/full_colinear_tab', quote = F, sep = "\t", row.names = F)
#write.table(LL_colinear_tab, 'data/colinear_LL_paralogs.tsv', quote = F, sep = "\t", row.names = F)
#write.table(AL_colinear_tab, 'data/colinear_AL_paralogs.tsv', quote = F, sep = "\t", row.names = F)
#write.table(XL_colinear_tab, 'data/colinear_XL_paralogs.tsv', quote = F, sep = "\t", row.names = F)



#+ scale_fill_manual(values=cbbPalette)
#cbbPalette <- c( "#009E73","#F0E442", "#E69F00","#CC79A7"  ,"#0072B2","#56B4E9", "#D55E00","#0072B2" , "#000000")
#ggplot(LL_colinear_tab, aes(y=mean_cov, x=anchored_scf, fill=block))  +
 # geom_point(position="dodge", stat="identity")+theme_classic()+theme(axis.text=element_text(size=14))

#boxplot(LL_colinear_tab$mean_cov~LL_colinear_tab$block )
#boxplot(LLb_colinear_tab$mean_cov~LLb_colinear_tab$block, add=T, col=c("blue") )

##trying to figure out how many A scaffolds are involved in L paralogs
ablocks<-AL_colinear_tab[AL_colinear_tab$chr=="A",]
ablocks<-ablocks[!duplicated(ablocks$anchored_scf), ]
##ALblocks come from 25 A scaffolds