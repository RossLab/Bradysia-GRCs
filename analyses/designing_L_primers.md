## Designing L primers

We would like to be able to study more carefully the two L chromosomes. For that it would be useful to have some markers on individual Ls. As a starting point we will take BUSCO genes that are X/A-L-L [tables/BUSCO_assigned.tsv](tables/BUSCO_assigned.tsv). We will then verify that the three scaffolds have coverages expected for the two different Ls.

```
# 2187at50557
2187at50557     NODE_25:144750-153487   2056    L       L-L-X   eIF-2-alpha kinase activator GCN1
2187at50557     NODE_376:45732-56314    2075    X       L-L-X   eIF-2-alpha kinase activator GCN1
2187at50557     NODE_611:39224-47254    1940    L       L-L-X   eIF-2-alpha kinase activator GCN1
```

```
NODE_25_length_173187_cov_11.987390: 142000 - 156000
NODE_611_length_60869_cov_15.293854: 43000 - 59000
NODE_376_length_71971_cov_106.929883: 37000 - 50000
```

~/generic_genomics/fasta2extract_by_header.py data/genome/genome.fasta NODE_25 > data/genome/genome_scf_NODE_25.fasta
~/generic_genomics/fasta2extract_by_header.py data/genome/genome.fasta NODE_611 > data/genome/genome_scf_NODE_611.fasta
~/generic_genomics/fasta2extract_by_header.py data/genome/genome.fasta NODE_376 > data/genome/genome_scf_NODE_376.fasta

~/generic_genomics/fasta2extract_by_nt_range.py data/genome/genome_scf_NODE_25.fasta 142000 156000 > data/genome/genome_scf_NODE_25_142000_156000.fasta
~/generic_genomics/fasta2extract_by_nt_range.py data/genome/genome_scf_NODE_611.fasta 43000 59000 > data/genome/genome_scf_NODE_611_43000_59000.fasta
~/generic_genomics/fasta2extract_by_nt_range.py data/genome/genome_scf_NODE_376.fasta 37000 50000 > data/genome/genome_scf_NODE_376_37000_50000.fasta


cat data/genome/genome_scf_NODE_376_37000_50000.fasta > data/genome/genome_2187at50557_sequences.fasta
cat data/genome/genome_scf_NODE_611_43000_59000.fasta >> data/genome/genome_2187at50557_sequences.fasta
cat data/genome/genome_scf_NODE_25_142000_156000.fasta >> data/genome/genome_2187at50557_sequences.fasta

### More genes

```R
BUSCO_genes <- read.table('tables/BUSCO_assigned.tsv', header = T, sep = '\t')
L_assignments <- read.table('data/L_scaffolds_chrom_assignment.tsv', header = T, sep = '\t')
scf_asn <- read.delim("data/scaffold_assignment_tab_full.tsv", header=T, stringsAsFactors = F)

row.names(L_assignments) <- L_assignments$scf
row.names(scf_asn) <- scf_asn$scf

sfc_names <- sapply(strsplit(BUSCO_genes$scf, ':'), function(x){ x[1] })

BUSCO_genes$L_asn <- L_assignments[sfc_names, 'L_asn']
# to save?? These few lines above will add L1 / L2

BUSCO_genes[BUSCO_genes$assignment_group == 'A-L-L', c('busco_id', 'scf', 'length', 'assignment', 'L_asn')]
# manually inspect
```

I got a list of target scaffolds out to `data/primers/list_of_scaffolds_to_extract`. Now I will get these scaffolds from the genome

```
~/generic_genomics/fasta2extract_by_list_of_headers.py data/genome/genome.fasta <(cat data/L-primers/list_of_scaffolds_to_extract | awk '{print $2}' | cut -f 1 -d ":") > data/L-primers/candidate_primer_sequences.fasta

```
