
# mapped kmers
sr_kmers <- read.table("data/Pacbio/6_kmermapping/table_of_mapped_kmers_spades.tsv", header = T, sep = '\t', stringsAsFactors = F)
sr_kmers$len <- sapply(strsplit(sr_kmers$id, '_'), function(x){ as.integer(x[4]) } )
sr_kmers$kmer_naive <- apply(sr_kmers[,c(2,3,4)], 1, function(x){ c('L', 'X', 'A')[which.max(x)] })
sr_kmers$score <- apply(sr_kmers[,c(2,3,4)], 1, function(x){ max(x) / sum(x) })

coverage_tab <- read.table('data/highcov_germ.tsv', stringsAsFactors = F)

# head(coverage_tab)
colnames(sr_kmers)[1] <- "name"

overlap_tab <- merge(coverage_tab, sr_kmers)

# the total 11 contigs coverage-wise L that are non-L by kmers
cov_yes_kmer_no <- overlap_tab[overlap_tab$kmer_naive == 'A',]
cov_yes_no_kmers_mapped <- coverage_tab[!coverage_tab$name %in% overlap_tab$name,]

# sum(c(cov_yes_kmer_no$len, cov_yes_no_kmers_mapped$len))
# 23775 nucleodites, not bad. Not bad at all. Given all this, it's quite clear that these will be sequences found both on L and A

# head(sr_kmers[sr_kmers$score < 0.9,])
# 'NODE_99_length_112107_cov_44.232839' is chimera

sr_kmers$AX <- sr_kmers$A + sr_kmers$X
sr_kmers$L_score <- sr_kmers$L / sr_kmers$len
sr_kmers$AX_score <- sr_kmers$AX / sr_kmers$len
# hist(sr_kmers$L_score)
# hist(sr_kmers$AX_score)

# library(hexbin)
# # Create hexbin object and plot
# h <- hexbin(data.frame(sr_kmers$L_score, sr_kmers$AX_score))
# plot(h)

L_kmer_candidates <- sr_kmers[sr_kmers$L_score > 0.8, c('name', 'L', 'len', 'L_score')]

# nrow(L_kmer_candidates)
# sum(L_kmer_candidates$len)

shared_candidates <- merge(coverage_tab, L_kmer_candidates)
kmer_specific_candidates <- L_kmer_candidates[!L_kmer_candidates$name %in% coverage_tab$name,]
coverage_specific_candidates <- coverage_tab[!coverage_tab$name %in% L_kmer_candidates$name,]

"both [Mbp]"
sum(shared_candidates$len) / 1e6
"kmers only [Mbp]"
sum(kmer_specific_candidates$len) / 1e6
"coverage only [Mbp]"
sum(coverage_specific_candidates$seqlen) / 1e6

shared_candidates <- shared_candidates[order(shared_candidates$len, decreasing = T),]
kmer_specific_candidates <- kmer_specific_candidates[order(kmer_specific_candidates$len, decreasing = T),]
coverage_specific_candidates <- coverage_specific_candidates[order(coverage_specific_candidates$seqlen, decreasing = T),]

write.table(shared_candidates[,'name'], file = "data/L-candidates-both.tsv",
            quote = F, sep = "\t", row.names = F, col.names = F)
write.table(kmer_specific_candidates[,'name'], file = "data/L-candidates-kmers-only.tsv",
            quote = F, sep = "\t", row.names = F, col.names = F)
write.table(coverage_specific_candidates[,'name'], file = "data/L-candidates-cov-only.tsv",
            quote = F, sep = "\t", row.names = F, col.names = F)

#### try to reconsider the normalization using number of UPPERCASE nucleotides
#
# scaffold_lengths <- read.table('data/genome/softmasking_per_scf.tsv', stringsAsFactors = F, header = T)
# scaffold_lengths[scaffold_lengths$name == 'NODE_60',]
#
# sr_kmers$name <- sapply(sr_kmers$name, function(scf){paste(unlist(strsplit(scf, "_"))[c(1,2)], collapse = "_")})
#
# srss_kmers = merge(sr_kmers, scaffold_lengths)
#
# sr_kmers$L_smscore <- sr_kmers$L / (srss_kmers$seqlen - srss_kmers$masked)
# sr_kmers$AX_smscore <- sr_kmers$AX / (srss_kmers$seqlen - srss_kmers$masked)
#
# h <- hexbin(data.frame(sr_kmers$L_smscore, sr_kmers$AX_smscore))
# plot(h)
#
# plot(sr_kmers$AX_smscore ~ sr_kmers$L_smscore, xlim = c(0, 3), ylim = c(0, 3))

# THIS IS A MESS
# because plenty of kmers are actually mapping on the softmasked positions and lead to nenormous mess up
