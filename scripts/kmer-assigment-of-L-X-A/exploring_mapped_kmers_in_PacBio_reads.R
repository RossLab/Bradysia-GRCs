mapped_kmers <- read.table("data/read_mpped_kmer.tsv", header = T, sep = '\t')

mapped_kmers$total <- rowSums(mapped_kmers[,c(2,3,4)])
mapped_kmers$max <- apply(mapped_kmers[,c(2,3,4)], 1, max)
mapped_kmers$score <- mapped_kmers$max / mapped_kmers$total
mapped_kmers$density <- mapped_kmers$total / mapped_kmers$len

hist(mapped_kmers$len, breaks = 60, main = 'histogram of read lengths')

hist(mapped_kmers$max, breaks = 120, main = 'histogram of total mapped kmers')

plot(mapped_kmers$total ~ mapped_kmers$len, pch = 20)

hist(mapped_kmers$score)
plot(mapped_kmers$score ~ mapped_kmers$len, pch = 20)
plot(mapped_kmers$score ~ mapped_kmers$total, pch = 20)


hist(, main = 'mapping kmer density', breaks = 120)

plot(mapped_kmers$total ~ mapped_kmers$len, pch = 20, xlim = c(0, 10000), ylim = c(0, 500))

plot(mapped_kmers$density ~ mapped_kmers$score, pch = 20)

# ~3416 out of 5000 :-/