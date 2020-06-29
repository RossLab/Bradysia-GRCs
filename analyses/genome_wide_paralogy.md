### Genome wide paralogy

The strategy is:
 - get all the transcripts
 - blast all vs all
 - extract all secondary hits

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

Now we need to do a magic. Generating a gtf file from the blast results (`data/genome_wide_paralogy/L_genes2PB_asm.blast`, `data/genome_wide_paralogy/AX_genes2ref_asm.blast`). Then reuse the all vs all protein blast and rerun MCScanX.

The gtf file needs to have 4 columns, scf, gene, from, to. I will use `blast_filter.py` from generic_genomics collection of scripts to get only the best hit per gene.

```
mkdir -p data/genome_wide_paralogy/anchored
ln -s ../proteins_all_vs_all.blast data/genome_wide_paralogy/anchored/Scop_anch_prot.blast
blast_filter.py data/genome_wide_paralogy/L_genes2PB_asm.blast 1 | awk '{ if( $3 > 95 ){ print $2 "\t" $1 "\t" $9 "\t" ($9 + $8 - $7) } }' > data/genome_wide_paralogy/anchored/Scop_anch_prot.gtf
blast_filter.py data/genome_wide_paralogy/AX_genes2ref_asm.blast 1 | awk '{ if( $3 > 95 ){ print $2 "\t" $1 "\t" $9 "\t" ($9 + $8 - $7) } }' >> data/genome_wide_paralogy/anchored/Scop_anch_prot.gtf
MCScanX -s 3 -m 100 -a data/genome_wide_paralogy/anchored/Scop_anch_prot
```



```
grep "^## " data/genome_wide_paralogy/anchored/Scop_anch_prot.collinearity | tr '=' ' ' | tr '&' ' ' | awk '{print $10 "\t" $11 "\t" $9 "\t" $5 "\t" $7}' > data/genome_wide_paralogy/Scop_anch_prot_collinearity.tsv
```

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

That's a bit crazy, there are genes that are in multiple collinear blocks that are completelly dominating the signal. These probably transposons that were annotated as genes, we should filter them out before we will carry on with the analysis.

Let's just generate a list of genes that are in collinear blocks, every gene will be there that many times as the number of block it is in

```bash
grep -v "^#" data/genome_wide_paralogy/anchored/Scop_anch_prot.collinearity | cut -f 2 > data/genome_wide_paralogy/colinear_genes.list
grep -v "^#" data/genome_wide_paralogy/anchored/Scop_anch_prot.collinearity | cut -f 3 >> data/genome_wide_paralogy/colinear_genes.list
```

and then in R, let's just make a list of those that are in more than 5 blocks (as the most obvious TE candidates).

```R
colinear_genes <- read.table('data/genome_wide_paralogy/colinear_genes.list', stringsAsFactors=F)$V1
in_number_of_blocks <- table(colinear_genes)
writeLines(names(in_number_of_blocks[in_number_of_blocks > 5]), 'data/genome_wide_paralogy/list_of_TE_candidates_to_filter.tsv')
```

Now we can use that list to filter out the blast and gff file to get filtered collinearity analysis and redo the analysis.

```
mkdir MCScanX data/genome_wide_paralogy/anchored_filtered
python3 scripts/genome_wide_paralogy/filter_gff_for_colinearity.py > data/genome_wide_paralogy/anchored_filtered/Scop_prot_anch_filt.gff
cp data/genome_wide_paralogy/anchored/*blast data/genome_wide_paralogy/anchored_filtered/Scop_prot_anch_filt.blast
MCScanX data/genome_wide_paralogy/anchored_filtered/Scop_prot_anch_filt
```

Now that's more like it. We get about 400 collinear blocks mostly A/X <-> L and L <-> L.

TODO
-> blast A/X to the PB assembly too (and do the analysis with PB fragmentation)
-> merge collinear blocks with the coverage table (do co-linear blocks have different coverages?)
-> quantify XA <-> XA, L <-> A and L <-> L (numbers and fractions)
-> verify that blasting made sense


##### Exonerate

This would be an alternative way how to "blast" genes. It does not run in parallel, so the way is to "divide and conquer" strategy, something like

```
exonerate --model protein2genome --target data/genome/ref/final_canu.fasta \
  --query data/genome/decomposed/genes_AX.fasta --showalignment no \
  --showcigar no --showvulgar no --bestn 1 --showquerygff no --showtargetgff yes \
  --targetchunkid $CHUNK --targetchunktotal 50 1> amphio_exonerate_chunk$CHUNK.out
```

However, I will try this only if blast won't work.

#### Based on scaffolds

`abacas.1.3.1.pl` does not really understand multifasta reference. Perhaps it's not the most optimal tool, but I will at least try to feed whole Illumina assembly to longest reference contig (`>contig_232`):

```
IL_ASM=data/genome/genome.fasta
abacas.1.3.1.pl -r data/genome/reference_genome_contig_232.fasta -q $IL_ASM -p nucmer
```

If that will produce some nice easily processable result, we might automatically separate the reference assembly into all scaffolds and do it. We might also map only the non-L scaffolds.