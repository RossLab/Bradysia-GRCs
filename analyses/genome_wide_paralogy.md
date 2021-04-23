### Genome wide paralogy

The strategy for paralog inference is:
 - get all the transcripts
 - blast all vs all
 - extract all reciprocal secondary hits
 - filter them by e value and alignment length
 - anochor the 

### getting the transcripts

```
gffread data/genome/annotation.gff3 -g data/genome/genome.fasta -w data/genome/transcripts.fasta -y data/genome/proteins.faa
```

Now convert transcripts into gene using `data/genome/transcripts2genes.map` (always taking the longest transcript as the gene sequence)

```
python3 scripts/genome_wide_paralogy/reduce_transcripts2genes.py data/genome/transcripts2genes.map data/genome/transcripts.fasta > data/genome/genes.fasta
python3 scripts/genome_wide_paralogy/reduce_transcripts2genes.py data/genome/transcripts2genes.map data/genome/proteins.faa -filter_list data/genome_wide_paralogy/list_of_transcripts_to_filter_out.list > data/genome/one_protein_per_gene.faa
```

now we got a `data/genome/genes.fasta` file we can finally selfblast.

### Genome wide paralogs using reciprocal blast

Using gene and protein sequences we infer the paralogy using reciprocal blast

**blasting**

```
GENES=data/genome/genes.fasta
makeblastdb -in $GENES -dbtype nucl
PROT=data/genome/one_protein_per_gene.faa
makeblastdb -in $PROT -dbtype prot

# blast all proteins vs all proteins
qsub -o logs -e logs -cwd -N selfblast -V -pe smp64 16 -b yes 'blastn -query '"$GENES"' -db '"$GENES"' -evalue 1e-10 -outfmt "6 std qlen slen" -num_threads 16 > data/genome_wide_paralogy/genes_all_vs_all.blast'
qsub -o logs -e logs -cwd -N selfblast -V -pe smp64 16 -b yes 'blastp -query '"$PROT"' -db '"$PROT"' -evalue 1e-10 -outfmt "6 std qlen slen" -num_threads 16 > data/genome_wide_paralogy/proteins_all_vs_all.blast'
```

**blast output processing**

Script `scripts/genome_wide_paralogy/reciprocal_blast.py` takes the blast output and generates a table of gene pairs with reciprocal hits (`_OG_pairs.tsv`) and  orthologous groups (`_OGs.tsv `). The table of pairs also contains % identity and length of alignment (taken from `gene1 gene2` blast record). Also parameter `-s` specifies the sequence similarity that is required to accept a blast alignment.

```
python3 scripts/genome_wide_paralogy/reciprocal_blast.py -s 60 data/genome_wide_paralogy/genes_all_vs_all.blast data/genome_wide_paralogy/nt_orthology
# takes a while (actually running it now and also take loads of memory, perhaps I need to do some optimization)
python3 scripts/genome_wide_paralogy/reciprocal_blast.py -s 60 data/genome_wide_paralogy/proteins_all_vs_all.blast data/genome_wide_paralogy/aa_orthology
```

The two files have also orthologous group IDs. These are just handles to easier work with the two files, but they don't correspond to each other between runs (if input or parameters change the orthologous group IDs as well).

### Nucleotide divergence of paralogs

We plot a quick and dirty histogram using this [this script](../scripts/genome_wide_paralogy/paralog_divergence_historgram.R):

![genome_wide_paralogy](https://user-images.githubusercontent.com/8181573/79567869-d6851e80-80ac-11ea-94bc-73aa3220a07b.png)

The divergence of L-L is (moreless) unimodel around 85% similarity, just like the L-A divergence that is around similar levels, suggesting an indipendent evolution of the two L copies shortly after the duplication.

### Collinearity analysis

```
# preparing annotation
awk '($3 == "mRNA") {
    OFS="\t";
    match($9, /ID=.+;Parent/);
    print $1, substr($9,RSTART+3,RLENGTH-10), $4, $5
}' data/genome/annotation.gff3 | sed 's/.t[1-9]//g' > data/genome_wide_paralogy/Scop_prot.gff

MCScanX -s 3 -m 100 -a data/genome_wide_paralogy/Scop_prot
```

reformating the output in something readable in R:

```
grep "^## " data/genome_wide_paralogy/Scop_prot.collinearity | tr '=' ' ' | tr '&' ' ' | awk '{print $10 "\t" $11 "\t" $9 "\t" $5 "\t" $7}' > data/genome_wide_paralogy/Scop_prot_collinearity.tsv
```

Ok, now getting in R the table

```
colinear_tab <- read.table('data/genome_wide_paralogy/Scop_prot_collinearity.tsv', stringsAsFactors = F, col.names = c('scf1', 'scf2', 'genes', 'score', 'eval'))
overlap_tab <- read.table('data/scaffold_assignment_tab_full.tsv', header = T, stringsAsFactors = F)

row.names(overlap_tab) <- overlap_tab$scf
colinear_tab$scf_1_ch <- overlap_tab[colinear_tab$scf1, 'assignments']
colinear_tab$scf_2_ch <- overlap_tab[colinear_tab$scf2, 'assignments']

paste(colinear_tab$scf_1_ch, colinear_tab$scf_2_ch)
```

There are plenty of LL blocks and several LA blocks. One LX block and no XA blocks. This is another evidence that we assembled the 2 L chromosomes, implying that the two L chomosomes are divergent from each other on nucleotide level but are indeed homologous. Note that although we detect relatively small number of synteny blocks (41 in total), Illumina assembly is very fragmented and hence there must be many more.

The next we try to redo the analysis using chromosome level asm of Sciara (maybe we can learn something from just remapping our non-L annotation on their assembly using exonerate, pasting together with our L and then rerunning `MCScanX`)

### anchoring Illumina genes

We are interested in gene evolution and the Illumina assembly has a substantially better BUSCO score than the PacBio assembly.

We can either try to get the annotation of the genome using something like [exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate), or directly blast the Illumina genes against PacBio (for L genes) and and the reference genome (for X and A). We could also try to do it on nucleotide level either directly by `nucmer` or `abacas` wrapper.

#### Based on genes

Extract genes that are AX-linked or L-linked. The AX-linked genes will be anchored on the reference assembly while the L-linked genes will be anchored on the PacBio assembly.

```
mkdir data/genome/decomposed
Rscript scripts/genome_wide_paralogy/generating_lists_of_genes.R
```

Getting the genes associated

```
python3 scripts/genome_wide_paralogy/reduce_transcripts2genes.py data/genome/transcripts2genes.map data/genome/transcripts.fasta -keep_list data/genome/decomposed/L_genes.list > data/genome/decomposed/genes_L.fasta
python3 scripts/genome_wide_paralogy/reduce_transcripts2genes.py data/genome/transcripts2genes.map data/genome/transcripts.fasta -keep_list data/genome/decomposed/X_genes.list > data/genome/decomposed/genes_X.fasta
python3 scripts/genome_wide_paralogy/reduce_transcripts2genes.py data/genome/transcripts2genes.map data/genome/transcripts.fasta -keep_list data/genome/decomposed/A_genes.list > data/genome/decomposed/genes_A.fasta
```

##### Blast

GRC genes to PacBio assmebly

```
L_GENES=data/genome/decomposed/genes_L.fasta
PB_asm=data/genome/pacbio_assembly.fasta
makeblastdb -in $PB_asm -dbtype nucl
qsub -o logs -e logs -cwd -N Lblast -V -pe smp64 16 -b yes "blastn -query $L_GENES -db $PB_asm -evalue 1e-10 -outfmt 6 -num_threads 16 > data/genome_wide_paralogy/L_genes2PB_asm.blast"
```

Non-L genes to the reference assembly

```
REF_asm=data/genome/ref/final_canu.fasta
mkdir -p data/genome/ref/
# the line bellow should be replaced by DL command once the asm will be online
ln -s /data/ross/mealybugs/analyses/pulling_sciara/Canu-assembly-Bcop_v1.0/Bcop_v1_assembly/final_canu.fasta $REF_asm
AX_GENES=data/genome/decomposed/genes_AX.fasta
cat data/genome/decomposed/genes_X.fasta data/genome/decomposed/genes_A.fasta > $AX_GENES

makeblastdb -in $REF_asm -dbtype nucl

qsub -o logs -e logs -cwd -N XAblast -V -pe smp64 16 -b yes "blastn -query $AX_GENES -db $REF_asm -evalue 1e-10 -outfmt 6 -num_threads 16 > data/genome_wide_paralogy/AX_genes2ref_asm.blast"
```

##### Generating blast-based annotation file

Now we need to do a magic. Generating a gff file from the blast results (`data/genome_wide_paralogy/L_genes2PB_asm.blast`, `data/genome_wide_paralogy/AX_genes2ref_asm.blast`). Then reuse the all vs all protein blast and rerun MCScanX.

The gff file needs to have 4 columns, scf, gene, from, to. I will use `blast_filter.py` from generic_genomics collection of scripts to get only the best hit per gene.

```
mkdir -p data/genome_wide_paralogy/anchored
ln -s ../proteins_all_vs_all.blast data/genome_wide_paralogy/anchored/Scop_anch_prot.blast
blast_filter.py data/genome_wide_paralogy/L_genes2PB_asm.blast 1 | awk '{ if( $3 > 95 ){ print $2 "\t" $1 "\t" $9 "\t" ($9 + $8 - $7) } }' > data/genome_wide_paralogy/anchored/Scop_anch_prot.gff
blast_filter.py data/genome_wide_paralogy/AX_genes2ref_asm.blast 1 | awk '{ if( $3 > 95 ){ print $2 "\t" $1 "\t" $9 "\t" ($9 + $8 - $7) } }' >> data/genome_wide_paralogy/anchored/Scop_anch_prot.gff
```

##### L orthologs mapped to the autosomes

This files contain 1. anchored our genes to the reference assembly and 2. chromosomal assignment of individual genes

```
data/genome_wide_paralogy/anchored/Scop_anch_prot.gff
data/nt_orthology_OG_pairs.tsv
```

We can merge them to figure out how many of the L genes map to anchored portion of the reference genome.

```{R}
orthologous_pairs <- read.table('data/nt_orthology_OG_pairs.tsv', header = T)

gene_anchoring <- read.table('data/genome_wide_paralogy/anchored/Scop_anch_prot.gff', header = F, col.names = c('scf', 'gene', 'from', 'to'))
rownames(gene_anchoring) <- gene_anchoring$gene

orthologous_pairs$anch1 <- gene_anchoring[orthologous_pairs$gene1, 'scf']
orthologous_pairs$anch2 <- gene_anchoring[orthologous_pairs$gene2, 'scf']

trans_asn <- read.delim("data/genome/gene.scaffold.map.tsv", header=F, stringsAsFactors = F, col.names = c('gene', 'scf'))
scf_asn <- read.delim("data/scaffold_assignment_tab_full.tsv", header=T, stringsAsFactors = F)

row.names(scf_asn) <- scf_asn$scf
trans_asn$assignments <- scf_asn[trans_asn$scf, 'assignments']
row.names(trans_asn) <- trans_asn$gene

orthologous_pairs[,c('chr1', 'chr2')] <- NA
orthologous_pairs$chr1 <- trans_asn[orthologous_pairs$gene1, 'assignments']
orthologous_pairs$chr2 <- trans_asn[orthologous_pairs$gene2, 'assignments']

AL_anchored <- orthologous_pairs[orthologous_pairs$chr1 == 'A' & orthologous_pairs$chr2 == 'L', 'anch1']
LA_anchored <- orthologous_pairs[orthologous_pairs$chr1 == 'L' & orthologous_pairs$chr2 == 'A', 'anch2']

sum(table(c(AL_anchored, LA_anchored))[c('contig_103')])
sum(table(c(AL_anchored, LA_anchored))[c('contig_171')])
sum(table(c(AL_anchored, LA_anchored))[c('contig_106', 'contig_144', 'contig_458', 'contig_5', 'contig_81', 'contig_84')])
```

We anchored 397 genes to chromosome II (`contig_103`), 80 to chromosome III (`contig_171`) and 374 gene to chromosome IV (contigs `contig_106,contig_144,contig_458,contig_5,contig_81,contig_84`). It does not seems that the L chromosome as whole is homologous to to any of the chromosomes, rather it represent bits and pieces from the whole genome. Furthermore, the the sizes of the anchored sequences to chromosomes II, III and IV (13.1, 5.4, and 46.4 Mbp), which is probably the reason why we have so "few" genes mapped to chromosome III. However, the L homologs / Mbp varies a lot between chromosomes (approx 30, 15 and 7 for the three chromosomes). We don't have any good explanation for that yet.

##### Collinearity analysis

```
MCScanX -s 3 -m 100 -a data/genome_wide_paralogy/anchored/Scop_anch_prot
```

and extract summary table from the full output table

```
grep "^## " data/genome_wide_paralogy/anchored/Scop_anch_prot.collinearity | tr '=' ' ' | tr '&' ' ' | awk '{print $10 "\t" $11 "\t" $9 "\t" $5 "\t" $7}' > data/genome_wide_paralogy/Scop_anch_prot_collinearity.tsv
```

now we explore the results in R

```{R}
Sciara_colinearity <- read.table('data/genome_wide_paralogy/Scop_anch_prot_collinearity.tsv', col.names = c('chr1', 'chr2', 'genes', 'score', 'eval'), stringsAsFactors = F )

length(unique(c(Sciara_colinearity$chr1, Sciara_colinearity$chr2)))
# 319
table(c(Sciara_colinearity$chr1, Sciara_colinearity$chr2))
# > PacBio assembly is still too fragmented to see big scale picture

head(Sciara_colinearity)
scfs_of_ch1_genes = unlist(apply(Sciara_colinearity, 1, function(x){ rep(x[1], as.numeric(x['genes'])) } ))
scfs_of_ch2_genes = unlist(apply(Sciara_colinearity, 1, function(x){ rep(x[2], as.numeric(x['genes'])) } ))

scaffolds_with_colinear_genes <- sort(table(c(scfs_of_ch1_genes, scfs_of_ch2_genes)))
scaffolds_with_colinear_genes <- data.frame(scf = names(scaffolds_with_colinear_genes),
                                            colinear_genes = as.vector(scaffolds_with_colinear_genes), stringsAsFactors = F)

# two_most <- Sciara_colinearity[Sciara_colinearity$chr1 %in% c('ctg171', 'ctg353') |
#                                Sciara_colinearity$chr2 %in% c('ctg171', 'ctg353')  ,]

hybrid_annotation <- read.table('data/genome_wide_paralogy/anchored/Scop_anch_prot.gff', stringsAsFactors = F, col.names = c('scf', 'gene', 'from', 'to'))
anchored_genes <- table(hybrid_annotation$scf)
anchored_genes <- data.frame(scf = names(anchored_genes), genes = as.vector(anchored_genes))
rownames(anchored_genes) <- anchored_genes$scf

scaffolds_with_colinear_genes$anchored_genes <- anchored_genes[scaffolds_with_colinear_genes$scf, 'genes']
```

##### Filtering using Christina's approach

Using approach in [this](https://github.com/RossLab/Sciara-L-chromosome/blob/master/scripts/paralog_divergence_prelim_analysis.R) script. I wrote a script that does this filtering directly on the blast file.

```
scripts/genome_wide_paralogy/filter_blast_for_colinearity.R
```

Now the filtered blast file (otherwise the same) is called `data/genome_wide_paralogy/genes_all_vs_all_filtered.blast`

```
mkdir data/genome_wide_paralogy/anchored_filtered
# python3 scripts/genome_wide_paralogy/filter_gff_for_colinearity.py > data/genome_wide_paralogy/anchored_filtered/Scop_prot_anch_filt.gff
ln -s `pwd`/data/genome_wide_paralogy/anchored/Scop_anch_prot.gff data/genome_wide_paralogy/anchored_filtered/Scop_anch_prot.gff
ln -s `pwd`/data/genome_wide_paralogy/proteins_all_vs_all_filtered.blast data/genome_wide_paralogy/anchored_filtered/Scop_anch_prot.blast

MCScanX data/genome_wide_paralogy/anchored_filtered/Scop_anch_prot
grep "^## " data/genome_wide_paralogy/anchored_filtered/Scop_anch_prot.collinearity | tr '=' ' ' | tr '&' ' ' | awk '{print $10 "\t" $11 "\t" $9 "\t" $5 "\t" $7}' > data/genome_wide_paralogy/Scop_anch_prot_filt.tsv
```

```{R}
collinear_table <- read.table('data/genome_wide_paralogy/Scop_anch_prot_filt.tsv', col.names = c('scf1', 'scf2', 'genes', 'score', 'eval'), stringsAsFactors = F)
X_scfs <- read.table('data/reference/X_contig_ID_john.txt', stringsAsFactors=F)[,1]

collinear_table[,c('chr1', 'chr2')] <- 'A'
collinear_table[grepl('ctg', collinear_table$scf1),'chr1'] <- 'L'
collinear_table[grepl('ctg', collinear_table$scf2),'chr2'] <- 'L'
collinear_table[collinear_table$scf1 %in% X_scfs,'chr1'] <- 'X'
collinear_table[collinear_table$scf2 %in% X_scfs,'chr2'] <- 'X'


table(paste(collinear_table$chr1, collinear_table$chr2))
# AX AX  AX L   L L
#     5    88    23

sum(collinear_table[collinear_table$chr1 == 'AX', 'genes'])
# 840
sum(collinear_table[collinear_table$chr1 == 'L', 'genes'])
# 172
```

Now that's more like it. We get 116 collinear blocks mostly A/X <-> L (88) carrying 860 genes and L <-> L (23) carrying 172 genes. Here is a script that will generate a table of genes with detailed information for downstream analses

```
python3 scripts/genome_wide_paralogy/merge_collinear_genes_with_coverages.py > data/genome_wide_paralogy/collinear_genes_full_table.tsv
```

Of 94 X <-> L scaffolds, 25 are anchored to individual chromosomes. Namely 12 to chromosome II (`contig_103`), 1 to chromosome III (`contig_171`) and 12 to chromosome IV (contigs `contig_106,contig_144,contig_458,contig_5,contig_81,contig_84`)

```
table(c(collinear_table$scf1, collinear_table$scf2))[c('contig_103', 'contig_171', 'contig_106', 'contig_144', 'contig_458', 'contig_5', 'contig_81', 'contig_84')]

# contig_103 contig_171 contig_106 contig_144 contig_458   contig_5       <NA>
#         15          1          4          2          5          2            
# contig_84
#         2

unique(c(collinear_table[collinear_table$scf1 %in% c('contig_103', 'contig_171', 'contig_106', 'contig_144', 'contig_458', 'contig_5', 'contig_81', 'contig_84'), 'scf2'],
  collinear_table[collinear_table$scf2 %in% c('contig_103', 'contig_171', 'contig_106', 'contig_144', 'contig_458', 'contig_5', 'contig_81', 'contig_84'), 'scf1']))

collinear_table[collinear_table$chr1 == 'X', ]

sum(collinear_table[collinear_table$scf1 %in% X_scfs, 'genes'])
# 119

all_links <- c(collinear_table$scf1, collinear_table$scf2)
multilink_Ls <- names(table(all_links[grepl('ctg', all_links)])[table(all_links[grepl('ctg', all_links)]) > 1])

collinear_table[collinear_table$scf1 %in% multilink_Ls | collinear_table$scf2 %in% multilink_Ls, ]
```

##### Evaluation of fragmentation effect

We mapped L to PB assembly and XA to the nearly chromosomal _Sciara_ reference, which is much more continuous than the PB assembly. Of course than detecting collinearity on A/X will be easier simply due to differences in fragmentation (N50). To get an idea about the effect of fragmentation we redo the analysis with PacBio assembly. We expect a dropout of blocks containing X/A genes.

```
PB_asm=data/genome/pacbio_assembly.fasta
AX_GENES=data/genome/decomposed/genes_AX.fasta
qsub -o logs -e logs -cwd -N AXblast -V -pe smp64 16 -b yes "blastn -query $AX_GENES -db $PB_asm -evalue 1e-10 -outfmt 6 -num_threads 16 > data/genome_wide_paralogy/AX_genes2PB_asm.blast"
```

TODO
 -> redo collinearity analysis
