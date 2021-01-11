library('ggplot2')

nwk <- read.table('tables/L-busco-phylogenies-summary.tsv', header = T)

nwk[,2][is.na(nwk[,2])] <- "unknown"
nwk[,3][is.na(nwk[,3])] <- "unknown"

head(nwk)

L_cecidomyiidae_branch_lengths <- c(nwk[nwk$GRC_1 == "cecidomyiidae", "GRC_1_len"], nwk[nwk$GRC_2 == "cecidomyiidae", "GRC_2_len"])
L_sciaridae_branch_lengths <- c(nwk[nwk$GRC_1 == "sciaridae", "GRC_1_len"], nwk[nwk$GRC_2 == "sciaridae", "GRC_2_len"])

teal <- rgb(76, 206, 175, maxColorValue = 255)
magenta <- rgb(206, 142, 218, maxColorValue = 255)

light_teal <- rgb(t(col2rgb(teal)), maxColorValue=255, alpha = 150)
light_magenta <- rgb(t(col2rgb(magenta)), maxColorValue=255, alpha = 150)

# hist(L_sciaridae_branch_lengths, breaks = 30, col = teal, probability = T, xlim = c(0,2))
# hist(L_cecidomyiidae_branch_lengths, breaks = 60, col = magenta, probability = T, add = T)

df <- data.frame(branch_length = c(L_sciaridae_branch_lengths, L_cecidomyiidae_branch_lengths),
 branch_type = c(rep('Sciaridae', length(L_sciaridae_branch_lengths)), rep('Cecidomyiidae', length(L_cecidomyiidae_branch_lengths))))

mu <- data.frame(branch_type = c('Sciaridae', 'Cecidomyiidae'), grp.mean = c(mean(L_sciaridae_branch_lengths), mean(L_cecidomyiidae_branch_lengths)))

p <- ggplot(df, aes(x=branch_length, color=branch_type, fill = branch_type)) +
  geom_density() +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=branch_type),
             linetype="dashed")

pdf('figures/GRC_genes_branch_length_distribution.pdf')

p + scale_color_manual(values = c(magenta, teal)) + scale_fill_manual(values=c(light_magenta, light_teal))

dev.off()

wilcox.test(L_cecidomyiidae_branch_lengths,L_sciaridae_branch_lengths)
