### Overview:
looking to see if GRC genes have a higher amino acid similarity to _M. destructor_ genes than genes in the core genome of _B. coprophila_. We want to look at this since most of the GRC BUSCO genes are in the Cecidomyiidae clade in phylogenies. But this analysis will take into account all the annotated genes.

### Analysis
all-by-all blast but with the Hessian fly annotation with the GRC genes appended. (then take out hessian fly only hits later).

I think the i5K assembly is from this paper" A BAC-based physical map of the Hessian fly genome anchored to polytene chromosomes".

hessian fly genome location: `/data/ross/mealybugs/analyses/Sciara-L-chromosome/data/m_destructor/genome`

file: `hf_OGS1.0.pep.fasta` this is protein file, also nt and gff

Number of genes: 22635

B. cop aa seq location:
`data/genome/one_protein_per_gene.faa`

Need to get a list of GRC geneID's and take those out of this file.

grc gene id's
`grc_gene_id.txt`

Taking grc genes and seq out of B. cop seq list
```
qsub -o logs -e logs -cwd -N grc_subset -V -pe smp64 1 -b yes 'seqtk subseq data/genome/one_protein_per_gene.faa data/grc_gene_id.txt > data/genome/grc_aa_geneseq.faa'
```
```
cat data/m_destructor/genome/hf_OGS1.0.pep.fasta data/genome/grc_aa_geneseq.faa > data/genome_wide_paralogy/grc_hf_aa_concat.fasta
```

Running blast and summarising
```
PROT=data/genome_wide_paralogy/grc_hf_aa_concat.fasta
makeblastdb -in $PROT -dbtype prot

# blast all proteins vs all proteins
qsub -o logs -e logs -cwd -N selfblast -V -pe smp64 16 -b yes 'blastp -query '"$PROT"' -db '"$PROT"' -evalue 1e-10 -outfmt "6 std qlen slen" -num_threads 16 > data/genome_wide_paralogy/hf.vs.grc_all_vs_all.blast'
```
```
qsub -o logs -e logs -cwd -N hf_blast -V -pe smp64 1 -b yes 'python3 scripts/genome_wide_paralogy/reciprocal_blast.py -s 40 data/genome_wide_paralogy/hf.vs.grc_all_vs_all.blast data/genome_wide_paralogy/hf.vs.grc_homology_40'
```

The rest of the analysis is in the script:
`scripts/hf.vs.grc_homologs.R`