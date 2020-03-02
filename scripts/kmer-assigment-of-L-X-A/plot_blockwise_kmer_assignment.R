
sr_kmers <- read.table("data/L-X-A-kmers/mapping/table_of_mapped_kmers_spades.tsv", header = T, sep = '\t', stringsAsFactors = F)
sr_kmers$len <- sapply(strsplit(sr_kmers$id, '_'), function(x){ as.integer(x[4]) } )

lr_kmer_blocks <- read.table('data/L-X-A-kmers/mapping/racon6pe_inspected_scfs_kmer_blocks.tsv', col.names = c('scf', 'from', 'to', 'L', 'X', 'A'), stringsAsFactors = F)
sr_kmer_blocks <- read.table('data/L-X-A-kmers/mapping/illumina_inspected_scfs_kmer_blocks.tsv', col.names = c('scf', 'from', 'to', 'L', 'X', 'A'), stringsAsFactors = F)

get_kmers_out <- function(x, type, gap){ c(sr_kmer_blocks[sr_kmer_blocks$scf == x, type], rep(NA, gap)) }
investigated_scfs = unique(sr_kmer_blocks$scf)

rownames(sr_kmers) <- sr_kmers$id
sr_kmers[investigated_scfs,]

per_scf_Ls <- lapply(investigated_scfs, get_kmers_out, 'L', 100)
per_scf_Xs <- lapply(investigated_scfs, get_kmers_out, 'X', 100)
per_scf_As <- lapply(investigated_scfs, get_kmers_out, 'A', 100)

plot(unlist(per_scf_Ls), col = 'yellow', ylim = c(0, 120))
points(unlist(per_scf_Xs), col = 'red')
points(unlist(per_scf_As), col = 'green')
for (border in cumsum(sapply(per_scf_Xs, length))){
    lines(c(border, border) - 50, c(0, 3))
}
legend('topright', bty = 'n', col = c('yellow', 'red', 'green'), pch = 20, c('L', 'X', 'A'))


# show those that have > 100 kmers in 100 long windows
sr_kmer_blocks[apply(sr_kmer_blocks[,4:6], 1, function(x){ any(x > 100) }), ]

get_lr_kmers_out <- function(x, type, gap){ c(lr_kmer_blocks[lr_kmer_blocks$scf == x, type], rep(NA, gap)) }
investigated_lr_scfs = unique(lr_kmer_blocks$scf)

per_scf_Ls <- lapply(investigated_lr_scfs, get_lr_kmers_out, 'L', 100)
per_scf_Xs <- lapply(investigated_lr_scfs, get_lr_kmers_out, 'X', 100)
per_scf_As <- lapply(investigated_lr_scfs, get_lr_kmers_out, 'A', 100)

plot(unlist(per_scf_Ls), col = 'yellow', ylim = c(0, 6000))
points(unlist(per_scf_Xs), col = 'red')
points(unlist(per_scf_As), col = 'green')
for (border in cumsum(sapply(per_scf_Xs, length))){
    lines(c(border, border) - 50, c(0, 5000))
}
legend('topright', bty = 'n', col = c('yellow', 'red', 'green'), pch = 20, c('L', 'X', 'A'))

hist(lr_kmer_blocks[lr_kmer_blocks$scf == "ctg33", 'A'], breaks = 60, col = 'green')
hist(lr_kmer_blocks[lr_kmer_blocks$scf == "ctg33", 'X'], breaks = 60, add = T, col = 'red')
hist(lr_kmer_blocks[lr_kmer_blocks$scf == "ctg33", 'L'], breaks = 20, add = T, col = 'yellow')

weird_illumina_kmer_blocks <- read.table('data/L-X-A-kmers/mapping/table_of_mapped_kmers_NODE_99_length_112107.tsv', col.names = c('scf', 'from', 'to', 'L', 'X', 'A'), stringsAsFactors = F)

get_kmers_out <- function(x, type, gap){ c(weird_illumina_kmer_blocks[weird_illumina_kmer_blocks$scf == x, type], rep(NA, gap)) }
investigated_scfs = unique(weird_illumina_kmer_blocks$scf)

per_scf_Ls <- lapply(investigated_scfs, get_kmers_out, 'L', 0)
per_scf_Xs <- lapply(investigated_scfs, get_kmers_out, 'X', 0)
per_scf_As <- lapply(investigated_scfs, get_kmers_out, 'A', 0)

png('figures/NODE_99_length_112107_chimera.png')
  plot(unlist(per_scf_Ls), col = 'yellow', ylim = c(0, 60), xlab = c('Position on scf [ 50bp ]'), ylab = 'Number of mapping kmers [kmers / 50bp window]', main = 'NODE_99_length_112107')
  points(unlist(per_scf_Xs), col = 'red')
  points(unlist(per_scf_As), col = 'green')
  legend('topright', bty = 'n', col = c('yellow', 'red', 'green'), pch = 20, c('L', 'X', 'A'))
dev.off()

hist(unlist(per_scf_Ls))
hist(unlist(per_scf_Xs))
