
# reading the BUSCO table
grcs <- read.table('tables/L-busco-phylogenies-summary.tsv', header = T)


grcs_L <- grcs[is.na(grcs$GRC_2), ]
single_GRC_busco <- table(grcs_L$GRC_1)
single_GRC_busco <- single_GRC_busco[c(1,3,2)]

# cecidomyiidae         other     sciaridae
#           247             6            86
 # 0.72861357 0.01769912 0.25368732

grcs_LL <- grcs[!is.na(grcs$GRC_2), ]
grcs_LL$both_grcs <- paste(grcs_LL$GRC_1, grcs_LL$GRC_2, sep = '-')

# these two are equivalent
grcs_LL$both_grcs[grcs_LL$both_grcs == 'sciaridae-cecidomyiidae'] <- 'cecidomyiidae-sciaridae'

double_GRC_busco <- table(grcs_LL$both_grcs)
# cecidomyiidae-cecidomyiidae         cecidomyiidae-other
#                          35                           2
#     cecidomyiidae-sciaridae         other-cecidomyiidae
#                          42                           1
#         sciaridae-sciaridae
#                           4

double_GRC_busco_other <- sum(double_GRC_busco[grepl('other', names(double_GRC_busco))])

double_GRC_busco <- c(double_GRC_busco[!grepl('other', names(double_GRC_busco))], double_GRC_busco_other)
names(double_GRC_busco)[4] <- 'other'

yl <- rgb(237,248,177, maxColorValue = 255)
gr <- rgb(127,205,187, maxColorValue = 255)
bl <- rgb(44,127,184, maxColorValue = 255)
pal <- c(bl, yl, 'grey', bl, gr, yl, 'grey')

pdf('figures/phylo_trees_topologies.pdf', width = 6, height = 4)

barplot(c(single_GRC_busco, double_GRC_busco), col = pal, space = c(0.1, 0.1, 0.1, 0.5, 0.1, 0.1, 0.1), cex.axis = 1, names.arg = F)

dev.off()
