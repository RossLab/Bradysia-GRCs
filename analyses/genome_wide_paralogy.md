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
python3 scripts/genome_wide_paralogy/reduce_transcripts2genes.py
```

now we got a `data/genome/genes.fasta` file we can finally selfblast.

### blasting

```
GENES=data/genome/genes.fasta
makeblastdb -in $GENES -dbtype nucl

# blast all proteins vs all proteins
qsub -o logs -e logs -cwd -N selfblast -V -pe smp64 16 -b yes "blastn -query $GENES -db $GENES -evalue 1e-10 -outfmt 6 -num_threads 16 | ~/generic_genomics/blast_filter.py 5 > data/genome_wide_paralogy/genes_all_vs_all.blast"
```