### Genome wide paralogy

The strategy is:
 - get all the transcripts (-M -K)
 - blast all vs all
 - extract all secondary hits

### getting the transcripts

```
gffread data/genome/annotation.gff3 -g data/genome/genome.fasta -w data/genome/transcripts.fasta -y data/genome/proteins.faa
```

Now convert transcripts into gene using `data/genome/transcripts2genes.map` (always taking the longest transcript as the gene sequence)

```
python3 scripts/genome_wide_paralogy/reducte_transcripts2genes.py
```

now we got a `data/genome/genes.fasta` file we can finally selfblast.

### blasting

