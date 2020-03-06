
# mapped kmers
sr_kmers <- read.table("data/L-X-A-kmers/mapping/table_of_mapped_kmers_spades.tsv", header = T, sep = '\t', stringsAsFactors = F)
colnames(sr_kmers)[1] <- "ID"

coverage_tab <- read.table('data/table.covdiff.germ.soma.txt', stringsAsFactors = F, header = T)
coverage_tab <- coverage_tab[coverage_tab$ID != "*",]
colnames(coverage_tab)[2] <- 'len'

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
overlap_tab$AX <- overlap_tab$A + overlap_tab$X
overlap_tab$L_score <- overlap_tab$L / overlap_tab$len
overlap_tab$AX_score <- overlap_tab$AX / overlap_tab$len

overlap_tab <- overlap_tab[order(overlap_tab$len, decreasing = T),]
write.table(overlap_tab, 'data/scaffold_assignment_tab_full.tsv', quote = F, sep = "\t", row.names = F)

L_kmer_candidates <- overlap_tab[overlap_tab$L_score > 0.8, c('ID', 'len', 'L_score', 'logdif2')]
coverage_candidates <- overlap_tab[overlap_tab$logdif2 > 2, c('ID', 'len', 'L_score', 'logdif2')]

shared_candidates <- merge(coverage_candidates, L_kmer_candidates)
kmer_specific_candidates <- L_kmer_candidates[!L_kmer_candidates$ID %in% coverage_candidates$ID,]
coverage_specific_candidates <- coverage_candidates[!coverage_candidates$ID %in% L_kmer_candidates$ID,]

"both [Mbp]"
sum(shared_candidates$len) / 1e6
"kmers only [Mbp]"
sum(kmer_specific_candidates$len) / 1e6
"coverage only [Mbp]"
sum(coverage_specific_candidates$len) / 1e6

shared_candidates <- shared_candidates[order(shared_candidates$len, decreasing = T),]
kmer_specific_candidates <- kmer_specific_candidates[order(kmer_specific_candidates$len, decreasing = T),]
coverage_specific_candidates <- coverage_specific_candidates[order(coverage_specific_candidates$len, decreasing = T),]

write.table(shared_candidates[,'ID'], file = "data/L-candidates-both.tsv",
            quote = F, sep = "\t", row.names = F, col.names = F)
write.table(kmer_specific_candidates[,'ID'], file = "data/L-candidates-kmers-only.tsv",
            quote = F, sep = "\t", row.names = F, col.names = F)
write.table(coverage_specific_candidates[,'ID'], file = "data/L-candidates-cov-only.tsv",
            quote = F, sep = "\t", row.names = F, col.names = F)
