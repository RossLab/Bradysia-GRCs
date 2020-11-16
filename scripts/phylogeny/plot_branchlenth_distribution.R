library('ggplot2')

nwk <- read.table('tables/L-busco-phylogenies-summary.tsv', header = T)

nwk[,2][is.na(nwk[,2])] <- "unknown"
nwk[,3][is.na(nwk[,3])] <- "unknown"

head(nwk)

L_cecidomyiidae_branch_lengths <- c(nwk[nwk$GRC_1 == "cecidomyiidae", "GRC_1_len"], nwk[nwk$GRC_2 == "cecidomyiidae", "GRC_2_len"])
L_sciaridae_branch_lengths <- c(nwk[nwk$GRC_1 == "sciaridae", "GRC_1_len"], nwk[nwk$GRC_2 == "sciaridae", "GRC_2_len"])

teal <- rgb(0.4098039, 0.8019608, 0.7431373)
magenta <- rgb(0.8411765, 0.7235294, 0.8803922)

# hist(L_sciaridae_branch_lengths, breaks = 30, col = teal, probability = T, xlim = c(0,2))
# hist(L_cecidomyiidae_branch_lengths, breaks = 60, col = magenta, probability = T, add = T)

df <- data.frame(branch_length = c(L_sciaridae_branch_lengths, L_cecidomyiidae_branch_lengths),
 branch_type = c(rep('Sciaridae', length(L_sciaridae_branch_lengths)), rep('Cecidomyiidae', length(L_cecidomyiidae_branch_lengths))))

mu <- data.frame(branch_type = c('Sciaridae', 'Cecidomyiidae'), grp.mean = c(mean(L_sciaridae_branch_lengths), mean(L_cecidomyiidae_branch_lengths)))

p <- ggplot(df, aes(x=branch_length, color=branch_type)) +
  geom_density() +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=branch_type),
             linetype="dashed")

p + scale_color_manual(values=c(magenta, teal))