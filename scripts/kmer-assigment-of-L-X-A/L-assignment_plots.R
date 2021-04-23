
overlap_tab <- read.table('data/scaffold_assignment_tab_full.tsv', header = T, stringsAsFactors = F)

overlap_tab$logdif2[is.infinite(overlap_tab$logdif2)] <- 20


# to Illumina spades assembly
sr_kmers <- read.table("data/L-X-A-kmers/mapping/table_of_mapped_kmers_spades.tsv", header = T, sep = '\t', stringsAsFactors = F)
sr_kmers$len <- sapply(strsplit(sr_kmers$id, '_'), function(x){ as.integer(x[4]) } )

# table_of_mapped_kmers.Lcandidates.tsv

expand_table <- function(mapped_kmers){
	mapped_kmers$total <- rowSums(mapped_kmers[,c(2,3,4)])
	mapped_kmers$max <- apply(mapped_kmers[,c(2,3,4)], 1, max)
	mapped_kmers$score <- mapped_kmers$max / mapped_kmers$len
	mapped_kmers$density <- mapped_kmers$total / mapped_kmers$len
	mapped_kmers$naive_assignment <- apply(mapped_kmers[,c(2,3,4)], 1, function(x){ c('L', 'X', 'A')[which.max(x)] })
	mapped_kmers
}

sr_kmers <- expand_table(sr_kmers)

add_one_ch <- function(mapped_kmers, ch = 'L', col = 'red', breaks = 120, feature = 'score'){
	hist(mapped_kmers[mapped_kmers$naive_assignment == ch, feature], breaks = breaks, add = T, col = col, freq = T)
}

pal <- c(rgb(229, 158, 37, maxColorValue = 255), rgb(32, 122, 95, maxColorValue = 255), rgb(108, 178, 216, maxColorValue = 255))

per_ch_hist <- function(mapped_kmers, feature, main = NA, xlim){
	hist(mapped_kmers[,feature], breaks = 120, freq = T, main = main, xlim = xlim, xlab = feature)
	add_one_ch(mapped_kmers, 'L', pal[1], feature = feature)
	add_one_ch(mapped_kmers, 'A', pal[2], feature = feature)
	add_one_ch(mapped_kmers, 'X', pal[3], feature = feature)
	legend('topleft', col = pal, pch = 20, c('GRC', 'A', 'X'), bty = 'n', cex = 1.4)
}

### Illumina Scores per naive assigments

png('figures/histograms_of_scores_per_group_IL.png')
	per_ch_hist(sr_kmers, 'score', 'Illumina max(L, X, A) / scf len', c(0, 1))
dev.off()

# png('figures/histograms_of_kmer_densities_per_group_IL.png')
# 	per_ch_hist(sr_kmers, 'density', 'Illumina sum(L, X, A) / scf len', c(0, 1))
# dev.off()

###

add_points_of_one_ch <- function(mapped_kmers, ch = 'L', col = 'red', feature1 = 'len', feature2 = 'score'){
	if ( feature1 == 'len') {
		f1 <- log10(mapped_kmers[mapped_kmers$naive_assignment == ch, 'len'])
	} else {
		f1 <- mapped_kmers[mapped_kmers$naive_assignment == ch, feature1]
	}

	f2 <- mapped_kmers[mapped_kmers$naive_assignment == ch, feature2]
	points(f1 ~ f2, pch = 20, col = col)
}

add.alpha <- function(cols, alpha) rgb(t(col2rgb(cols)/255), alpha = alpha)
pal <- add.alpha(pal, 0.3)

plot_dot_plot <- function(mapped_kmers, feature1, feature2){
	if ( feature1 == 'len') {
		f1 <- log10(mapped_kmers[, 'len'])
	} else {
		f1 <- mapped_kmers[, feature1]
	}
	plot(f1 ~ mapped_kmers[,feature2], pch = 20, cex = 1.2, ylab = feature1, xlab = feature2)
	add_points_of_one_ch(mapped_kmers, 'L', pal[1], feature1 = feature1, feature2 = feature2)
	add_points_of_one_ch(mapped_kmers, 'A', pal[2], feature1 = feature1, feature2 = feature2)
	add_points_of_one_ch(mapped_kmers, 'X', pal[3], feature1 = feature1, feature2 = feature2)
	legend('topleft', col = pal, pch = 20, c('L', 'A', 'X'), bty = 'n')
}

png('figures/length_vs_score_IL.png')
	plot_dot_plot(sr_kmers, 'len', 'score')
# density and score are correlated
dev.off()
