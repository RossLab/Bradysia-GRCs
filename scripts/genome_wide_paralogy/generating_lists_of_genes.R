trans_asn <- read.delim("data/genome/gene.scaffold.map.tsv", header=F, stringsAsFactors = F, col.names = c('gene', 'scf'))
scf_asn <- read.delim("data/scaffold_assignment_tab_full.tsv", header=T, stringsAsFactors = F)

row.names(scf_asn) <- scf_asn$scf
trans_asn$assignments <- scf_asn[trans_asn$scf, 'assignments']
# removing NA values
trans_asn <- trans_asn[!is.na(trans_asn$assignments),]

# write lists
for (asn in c('L', 'X', 'A')){
    write(trans_asn[trans_asn$assignments == asn,'gene'],
          paste0("data/genome/decomposed/", asn ,"_genes.list"))
}

for (asn in c('L', 'X', 'A')){
    write(scf_asn[scf_asn$assignments == asn,'scf'],
          paste0("data/genome/decomposed/", asn ,"_scfs.list"))
}