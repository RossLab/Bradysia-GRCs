#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(hexbin)
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

kmer_dump <- read.table(args[1], col.names = c('heads', 'testes'))

png('figures/sciara_head_testes_2d_hist.png')
	hexbinplot(testes ~ heads, data=kmer_dump, colramp=rf)
dev.off()