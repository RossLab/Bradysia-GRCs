
buscos <- read.table('tables/L-busco-phylogenies-summary-L-in-Sciaridae_only.tsv', header = T)

assignments <- strsplit(buscos[, 2], ',')
types <- strsplit(buscos[, 3], ',')

just_Ls <- lapply(1:nrow(buscos), function(i){types[[i]][assignments[[i]] == 'L']})

table(unlist(just_Ls))