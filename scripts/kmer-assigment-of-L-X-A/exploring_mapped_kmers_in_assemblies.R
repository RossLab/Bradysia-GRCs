


# to Illumina spades assembly
sr_kmers <- read.table("data/Pacbio/6_kmermapping/table_of_mapped_kmers_spades.tsv", header = T, sep = '\t', stringsAsFactors = F)
sr_kmers$len <- sapply(strsplit(sr_kmers$id, '_'), function(x){ as.integer(x[4]) } )
# to PacBio assembly
lr_kmers <- read.table("data/Pacbio/6_kmermapping/table_of_mapped_kmers_PacBio.tsv", header = T, sep = '\t', stringsAsFactors = F)
lr_contig_lengths <- read.table('data/Pacbio/6_kmermapping/racon6pe_contig_lengths.tsv', sep = '\t', header = F, stringsAsFactors = F, col.names = c('id', 'len'))
lr_kmers <- merge(lr_kmers, lr_contig_lengths)
# sorting them by scaffold name
lr_kmers <- lr_kmers[order(lr_kmers$len, decreasing = T),]

expand_table <- function(mapped_kmers){
	mapped_kmers$total <- rowSums(mapped_kmers[,c(2,3,4)])
	mapped_kmers$max <- apply(mapped_kmers[,c(2,3,4)], 1, max)
	mapped_kmers$score <- mapped_kmers$max / mapped_kmers$total
	mapped_kmers$density <- mapped_kmers$total / mapped_kmers$len
	mapped_kmers$naive_assignment <- apply(mapped_kmers[,c(2,3,4)], 1, function(x){ c('L', 'X', 'A')[which.max(x)] })
	mapped_kmers
}

sr_kmers <- expand_table(sr_kmers)
lr_kmers <- expand_table(lr_kmers)

####

pal <- c(rgb(0.9,0.2,0.2, alpha = 0.5), rgb(0.2,0.2,0.9, alpha = 0.5))

png('figures/histograms_of_scores.png')

hist(lr_kmers$score, breaks = 120, main = 'max(L, X, A) / sum(L, X, A)', freq = F, col = pal[2], xlim = c(0.8,1), xlab = c('Score'))
hist(sr_kmers$score, breaks = 120, freq = F, add = T, col = pal[1])
legend('topleft', bty = 'n', col = pal, c('Illumina', 'PacBio'), pch = 20, cex = 1.4)

dev.off()

### PacBio Scores per naive assigments

add_one_ch <- function(mapped_kmers, ch = 'L', col = 'red', breaks = 120, feature = 'score'){
	hist(mapped_kmers[mapped_kmers$naive_assignment == ch, feature], breaks = breaks, add = T, col = col, freq = T)
}

per_ch_hist <- function(mapped_kmers, feature, main = NA, xlim){
	hist(mapped_kmers[,feature], breaks = 120, freq = T, main = main, xlim = xlim, xlab = feature)
	add_one_ch(mapped_kmers, 'L', 'yellow', feature = feature)
	add_one_ch(mapped_kmers, 'A', 'green', feature = feature)
	add_one_ch(mapped_kmers, 'X', 'red', feature = feature)
	legend('topleft', col =  c('yellow', 'green', 'red'), pch = 20, c('L', 'A', 'X'), bty = 'n', cex = 1.4)
}



png('figures/histograms_of_scores_per_group_PB.png')
	per_ch_hist(lr_kmers, 'score', 'PacBio max(L, X, A) / sum(L, X, A)', c(0.8, 1))
dev.off()


### Illumina Scores per naive assigments

png('figures/histograms_of_scores_per_group_IL.png')
	per_ch_hist(sr_kmers, 'score', 'Illumina max(L, X, A) / sum(L, X, A)', c(0.8, 1))
dev.off()

###

png('figures/histograms_of_kmer_densities_per_group_IL.png')
	per_ch_hist(sr_kmers, 'density', 'Illumina sum(L, X, A) / scf len', c(0, 1))
dev.off()

###

png('figures/histograms_of_kmer_densities_per_group_PB.png')
	per_ch_hist(lr_kmers, 'density', 'PacBIo sum(L, X, A) / scf len', c(0, 1))
# wow, the X/A scfs are a mess -> so low densities. Is it because of 
dev.off()

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

plot_dot_plot <- function(mapped_kmers, feature1, feature2){
	if ( feature1 == 'len') {
		f1 <- log10(mapped_kmers[, 'len'])		
	} else {
		f1 <- mapped_kmers[, feature1]		
	}
	plot(f1 ~ mapped_kmers[,feature2], pch = 20, cex = 1.2, ylab = feature1, xlab = feature2)
	add_points_of_one_ch(mapped_kmers, 'L', 'yellow', feature1 = feature1, feature2 = feature2)
	add_points_of_one_ch(mapped_kmers, 'A', 'green', feature1 = feature1, feature2 = feature2)
	add_points_of_one_ch(mapped_kmers, 'X', 'red', feature1 = feature1, feature2 = feature2)
	legend('topleft', col =  c('yellow', 'green', 'red'), pch = 20, c('L', 'A', 'X'), bty = 'n')
}

png('figures/length_vs_score_IL.png')
	plot_dot_plot(sr_kmers, 'len', 'score')
# density and score are correlated
dev.off()

png('figures/density_vs_score_IL.png')
	plot_dot_plot(sr_kmers, 'density', 'score')
# density and score are correlated
dev.off()

png('figures/length_vs_density_IL.png')
	plot_dot_plot(sr_kmers, 'len', 'density')
dev.off()

png('figures/length_vs_score_PB.png')
	plot_dot_plot(lr_kmers, 'len', 'score')
# density and score are correlated
dev.off()

png('figures/density_vs_score_PB.png')
	plot_dot_plot(lr_kmers, 'density', 'score')
# density and score are correlated
dev.off()

png('figures/length_vs_density_PB.png')
	plot_dot_plot(lr_kmers, 'len', 'density')
# also true for PB, but a bit less strongly
dev.off()


sapply(c('X', 'L', 'A'), function(x){ sum(sr_kmers[sr_kmers$naive_assignment == x, 'len']) } )
sapply(c('X', 'L', 'A'), function(x){ sum(sr_kmers[lr_kmers$naive_assignment == x, 'len']) } )