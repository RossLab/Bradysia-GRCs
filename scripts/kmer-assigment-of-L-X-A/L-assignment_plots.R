
library(hexbin)

overlap_tab <- read.table('data/scaffold_assignment_tab_full.tsv', header = T, stringsAsFactors = F)

overlap_tab$logdif2[is.infinite(overlap_tab$logdif2)] <- 20

# # Create hexbin object and plot
h <- hexbin(data.frame(overlap_tab$L_score, overlap_tab$logdif2))

png('figures/coverage_kmers_density.png')
  plot(h, xlab = "L score", ylab = "log2 coverage difference", main = 'density of coverage and kmer method')
dev.off()

# log of the densities. The plot handles only integers, so I add 1 to all non-zero values to be at least 1 on the plot (otherwise log2(1) would round to 0)
h@count <- log2(h@count + (h@count > 0))
png('figures/coverage_kmers_density_log.png')
  plot(h, xlab = "L score", ylab = "log2 coverage difference", main = 'log density of coverage and kmer method')
dev.off()

# TESTING THE COVERAGE RATIO IDEA, it's kind of messy and it's hard to incorporate the density there explicitly
# overlap_tab$log2diff_kmer_ratio <- log2(overlap_tab$L / (overlap_tab$A + overlap_tab$X) * (overlap_tab$L_score + overlap_tab$AX_score))
# overlap_tab$log2diff_kmer_ratio[is.infinite(overlap_tab$log2diff_kmer_ratio)] <- 20
# overlap_tab$log2diff_kmer_ratio_normalized <- overlap_tab$log2diff_kmer_ratio
#
# h <- hexbin(data.frame(overlap_tab$log2diff_kmer_ratio_normalized, overlap_tab$logdif2))
# plot(h)
#
# h@count <- log2(h@count + (h@count > 0))
# plot(h)
#
# h <- hexbin(data.frame(overlap_tab$log2diff_kmer_ratio * overlap_tab$density, overlap_tab$L_score))
# plot(h)

# h@count <- log2(h@count + (h@count > 0))
# plot(h)

# hist(abs(log2(overlap_tab$L_score)) * overlap_tab$logdif2, breaks = 16000, xlim = c(-5, 5))


#### try to reconsider the normalization using number of UPPERCASE nucleotides
#
# scaffold_lengths <- read.table('data/genome/softmasking_per_scf.tsv', stringsAsFactors = F, header = T)
# scaffold_lengths[scaffold_lengths$name == 'NODE_60',]
#
# overlap_tab$name <- sapply(overlap_tab$name, function(scf){paste(unlist(strsplit(scf, "_"))[c(1,2)], collapse = "_")})
#
# srss_kmers = merge(overlap_tab, scaffold_lengths)
#
# overlap_tab$L_smscore <- overlap_tab$L / (srss_kmers$seqlen - srss_kmers$masked)
# overlap_tab$AX_smscore <- overlap_tab$AX / (srss_kmers$seqlen - srss_kmers$masked)
#
# h <- hexbin(data.frame(overlap_tab$L_smscore, overlap_tab$AX_smscore))
# plot(h)
#
# plot(overlap_tab$AX_smscore ~ overlap_tab$L_smscore, xlim = c(0, 3), ylim = c(0, 3))

# THIS IS A MESS
# because plenty of kmers are actually mapping on the softmasked positions and lead to nenormous mess up
