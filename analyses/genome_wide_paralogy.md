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

### blasting

```
GENES=data/genome/genes.fasta
makeblastdb -in $GENES -dbtype nucl
PROT=data/genome/one_protein_per_gene.faa
makeblastdb -in $PROT -dbtype prot

# blast all proteins vs all proteins
qsub -o logs -e logs -cwd -N selfblast -V -pe smp64 16 -b yes "blastn -query $GENES -db $GENES -evalue 1e-10 -outfmt 6 -num_threads 16 | ~/generic_genomics/blast_filter.py 5 > data/genome_wide_paralogy/Scop_prot.blast"
qsub -o logs -e logs -cwd -N selfblast -V -pe smp64 16 -b yes "blastp -query $PROT -db $PROT -evalue 1e-10 -outfmt 6 -num_threads 16 | ~/generic_genomics/blast_filter.py 5 > data/genome_wide_paralogy/proteins_all_vs_all.blast"
```

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

There are two things we could do more:
 - adjust this analysis for paralogs between L scaffolds and figure how many synteny breaks are there (how many paralogs are in different orders rather than the same). This will give us an idea about preservation of the order of genes on the two Ls
 - alternatively we could try to redo the analysis using chromosome level asm of Sciara (maybe we can learn something from just remapping our non-L annotation on their assembly using exonerate, pasting together with our L and then rerunning `MCScanX`) 