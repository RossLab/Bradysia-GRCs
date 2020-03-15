
scf_assignment_tab <- read.table('data/scaffold_assignment_tab_full.tsv', stringsAsFactors = F, header = T)
row.names(scf_assignment_tab) <- as.vector(sapply(scf_assignment_tab$ID, function(x) { x = unlist(strsplit(x, "_")); paste0(x[1], '_', x[2]) } ))

BUSCO_full_tab <- read.table('data/busco.sciara.masked.spades.long_full_table.tsv',
                             stringsAsFactors = F, header = F, sep = '\t', quote="\"", fill = T,
                             col.names = c('busco_id','status', 'scf', 'start', 'end', 'score', 'length', 'URL', 'fction'))


# filtering out 12 missing genes
BUSCO_full_tab <- BUSCO_full_tab[BUSCO_full_tab$status != "Missing",]

gather_duplicates <- function(gene){
  scaffolds <- BUSCO_full_tab[BUSCO_full_tab$busco_id == gene,'scf']
  sapply(strsplit(scaffolds, ":"), function(x){ x[1]} )
}

# getting BUSCO -> group (i.e. chromosomal topology)
BUSCO_genes <- unique(BUSCO_full_tab$busco_id)
list_of_paralogs <- lapply(BUSCO_genes, gather_duplicates)
busco_assignment_groups <- sapply(list_of_paralogs, function(x){paste(sort(scf_assignment_tab[x, 'assignments']), collapse = "-")} )
names(busco_assignment_groups) <- BUSCO_genes
# BUSCO_genes[busco_assignment_groups == 'AX-AX-L-L-L']
# BUSCO_full_tab[BUSCO_full_tab$busco_id == '57373at50557',]

scfs <- sapply(strsplit(BUSCO_full_tab$scf, ":"), function(x){ x[1]} )
BUSCO_full_tab$assignment <- scf_assignment_tab[scfs,'assignments']
BUSCO_full_tab$assignment_group <- busco_assignment_groups[BUSCO_full_tab$busco_id]

png('figures/BUSCO_scores.png')
  hist(BUSCO_full_tab[BUSCO_full_tab$assignment == 'AX', 'score'], breaks = 60, col = 'purple', main = "BUSCO score vs chromosomal assignment", xlab = 'BUSCO score')
  hist(BUSCO_full_tab[BUSCO_full_tab$assignment == 'L', 'score'], breaks = 30, col = 'yellow', add = T)
  legend('topright', bty = 'n', pch = 20, c('GRC', 'AX'), col = c('yellow', 'purple'))
dev.off()

output_table <- BUSCO_full_tab[,c('busco_id','scf','length', 'assignment', 'assignment_group','fction')]
write.table(output_table, 'tables/BUSCO_assigned.tsv', quote = F, sep = "\t", row.names = F)
# To load this tab to R
# read.table('tables/BUSCO_assigned.tsv', header = T, sep = '\t', stringsAsFactors = F, quote = "\"")
