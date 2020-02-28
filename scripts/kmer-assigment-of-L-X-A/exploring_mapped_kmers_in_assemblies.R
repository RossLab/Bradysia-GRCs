


# to Illumina spades assembly
sr_kmers <- read.table("data/Pacbio/6_kmermapping/table_of_mapped_kmers_spades.tsv", header = T, sep = '\t', stringsAsFactors = F)
sr_kmers$len <- sapply(strsplit(sr_kmers$id, '_'), function(x){ as.integer(x[4]) } )

# to PacBio assembly
lr_kmers <- read.table("data/Pacbio/6_kmermapping/table_of_mapped_kmers_PacBio.tsv", header = T, sep = '\t', stringsAsFactors = F)
lr_contig_lengths <- read.table('data/Pacbio/6_kmermapping/racon6pe_contig_lengths.tsv', sep = '\t', header = F, stringsAsFactors = F, col.names = c('id', 'len'))
lr_kmers <- merge(lr_kmers, lr_contig_lengths)
# sorting them by scaffold name
lr_kmers <- lr_kmers[order(lr_kmers$len, decreasing = T),]

# to John's assembly
j_kmers <- read.table("data/Pacbio/6_kmermapping/table_of_mapped_kmers.johnassembly.tsv", header = T, sep = '\t', stringsAsFactors = F)
j_contig_lengths <- read.table('data/Pacbio/6_kmermapping/johnassembly_contig_lengths.tsv', sep = '\t', header = F, stringsAsFactors = F, col.names = c('id', 'len'))
j_kmers <- merge(j_kmers, j_contig_lengths)
# sorting them by scaffold name
j_kmers <- j_kmers[order(j_kmers$len, decreasing = T),]

jlc_kmers <- read.table("data/Pacbio/6_kmermapping/table_of_mapped_kmers.Lcandidates.tsv", header = T, sep = '\t', stringsAsFactors = F)
jlc_contig_lengths <- read.table('data/Pacbio/6_kmermapping/Lcandidate_contig_lengths.tsv', sep = '\t', header = F, stringsAsFactors = F, col.names = c('id', 'len'))
jlc_kmers <- merge(jlc_kmers, jlc_contig_lengths)
# sorting them by scaffold name
jlc_kmers <- jlc_kmers[order(jlc_kmers$len, decreasing = T),]

# table_of_mapped_kmers.Lcandidates.tsv

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
j_kmers <- expand_table(j_kmers)
jlc_kmers <- expand_table(jlc_kmers)

####

pal <- c(rgb(0.9,0.2,0.2, alpha = 0.5), rgb(0.2,0.2,0.9, alpha = 0.5), rgb(0.9,0.2,0.9, alpha = 0.5), rgb(0.2,0.9,0.5, alpha = 0.5))

png('figures/histograms_of_scores.png')

hist(lr_kmers$score, breaks = 120, main = 'max(L, X, A) / sum(L, X, A)', freq = F, col = pal[2], xlim = c(0.8,1), xlab = c('Score'))
hist(sr_kmers$score, breaks = 120, freq = F, add = T, col = pal[1])
hist(j_kmers$score, breaks = 120, freq = F, add = T, col = pal[3])
hist(jlc_kmers$score, breaks = 120, freq = F, add = T, col = pal[4])
legend('topleft', bty = 'n', col = pal, c('Illumina', 'PacBio', "John asm", "L candidates"), pch = 20, cex = 1.4)

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

png('figures/histograms_of_scores_per_group_J.png')
	per_ch_hist(j_kmers, 'score', 'John max(L, X, A) / sum(L, X, A)', c(0.8, 1))
dev.off()

png('figures/histograms_of_scores_per_group_JLC.png')
	per_ch_hist(jlc_kmers, 'score', 'John max(L, X, A) / sum(L, X, A)', c(0.8, 1))
dev.off()
###

png('figures/histograms_of_kmer_densities_per_group_IL.png')
	per_ch_hist(sr_kmers, 'density', 'Illumina sum(L, X, A) / scf len', c(0, 1))
dev.off()

# sum(jlc_kmers[jlc_kmers$naive_assignment == "L", 'len'])
# [1] 3303504
# sum(jlc_kmers[jlc_kmers$naive_assignment != "L", 'len'])
# [1] 9412304

###

png('figures/histograms_of_kmer_densities_per_group_PB.png')
	per_ch_hist(lr_kmers, 'density', 'PacBIo sum(L, X, A) / scf len', c(0, 1))
# wow, the X/A scfs are a mess -> so low densities. Is it because of
dev.off()

png('figures/histograms_of_kmer_densities_per_group_J.png')
	per_ch_hist(j_kmers, 'density', 'PacBIo sum(L, X, A) / scf len', c(0, 1))
dev.off()

png('figures/histograms_of_kmer_densities_per_group_J_Lcandidates.png')
	per_ch_hist(jlc_kmers, 'density', 'L candidates sum(L, X, A) / scf len', c(0, 1))
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


plot_dot_plot(j_kmers, 'len', 'score')
plot_dot_plot(j_kmers, 'density', 'score')
plot_dot_plot(j_kmers, 'len', 'density')
# also true for PB, but a bit less strongly

sapply(c('X', 'L', 'A'), function(x){ round(sum(sr_kmers[sr_kmers$naive_assignment == x, 'len']) / 1e6, 2) } )
sapply(c('X', 'L', 'A'), function(x){ round(sum(lr_kmers[lr_kmers$naive_assignment == x, 'len']) / 1e6, 2) } )
sapply(c('X', 'L', 'A'), function(x){ round(sum(j_kmers[j_kmers$naive_assignment == x, 'len']) / 1e6, 2) } )

head(j_kmers[,2:4])

# expand_get_combinations <- function(mapped_kmers){
# 	mapped_kmers$second_most_frequent <- apply(mapped_kmers[,c(2,3,4)], 1, function(x){ c('L', 'X', 'A')[-c(which.min(x), which.max(x))] })
# 	mapped_kmers$fist_second_diff <- apply(mapped_kmers[,c(2,3,4)], 1, function(x){ x_sorted <- sort(x, decreasing=T); (x_sorted[1] - x_sorted[2]) })
# 	mapped_kmers
# }
#
# sr_kmers <- expand_get_combinations(sr_kmers)
# lr_kmers <- expand_get_combinations(lr_kmers)
# j_kmers <- expand_get_combinations(j_kmers)
#
# sr_kmers$pair <- paste(sr_kmers$naive_assignment, sr_kmers$second_most_frequent)
# lr_kmers$pair <- paste(lr_kmers$naive_assignment, lr_kmers$second_most_frequent)
# j_kmers$pair <- paste(j_kmers$naive_assignment, j_kmers$second_most_frequent)
#
# boxplot(j_kmers$score ~ j_kmers$pair)
# boxplot(lr_kmers$score ~ lr_kmers$pair)
# boxplot(sr_kmers$score ~ sr_kmers$pair)
#
# pairs <- c('A L', 'A X', 'L A', 'L X', 'X A', 'X L')
# sapply(pairs, function(x) { sum(sr_kmers[sr_kmers$pair == x, 'len']) }) / 1e6
# sapply(pairs, function(x) { sum(lr_kmers[lr_kmers$pair == x, 'len']) }) / 1e6
# sapply(pairs, function(x) { sum(j_kmers[j_kmers$pair == x, 'len']) }) / 1e6
