input_blast_file <- 'data/genome_wide_paralogy/proteins_all_vs_all.blast'
output_blast_file <- 'data/genome_wide_paralogy/proteins_all_vs_all_filtered.blast'

blast_frame <- read.table(input_blast_file, sep = '\t', stringsAsFactors = F)
colnames(blast_frame) <- c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen')

# filter_paralog_cov_tab<-paralog_cov_tab[paralog_cov_tab$lendif1<=2,]
# filter_paralog_cov_tab<-filter_paralog_cov_tab[filter_paralog_cov_tab$lendif2<=2,]
two_fold_size_filter <- (blast_frame$qlen / blast_frame$slen) <= 2 & (blast_frame$slen / blast_frame$qlen) <= 2
# filter_paralog_cov_tab<-filter_paralog_cov_tab[filter_paralog_cov_tab$percentaln1 >0.7 , ]
# filter_paralog_cov_tab<-filter_paralog_cov_tab[filter_paralog_cov_tab$percentaln2 >0.7 , ]
alignment_cov_filt <- (blast_frame$length / blast_frame$slen) > 0.7 & (blast_frame$length / blast_frame$qlen) > 0.7

blast_frame <- blast_frame[two_fold_size_filter,]
write.table(blast_frame, output_blast_file, sep = '\t', quote = F, col.names = F, row.names = F)
