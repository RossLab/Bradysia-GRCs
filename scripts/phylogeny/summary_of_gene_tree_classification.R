
buscos <- read.table('tables/classification_of_GRC_genes_places_in_sciaridae.tsv', header = T)

assignments <- strsplit(buscos[, 2], ',')
types <- strsplit(buscos[, 3], ',')

just_Ls <- lapply(1:nrow(buscos), function(i){types[[i]][assignments[[i]] == 'L']})

table(unlist(just_Ls))