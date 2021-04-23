### Estimation of GRC size using kmer spectra

Is it possible to get the L size using kmers?

```
scripts/kmer-assigment-of-L-X-A/dump2hist.py data/L-X-A-kmers/kmer_db/merged_k27.dump data/L-X-A-kmers/L_kmers_k27.hist
```

generates the histogram that I will feed to genomescope

```
genomescope.R -i data/L-X-A-kmers/L_kmers_k27.hist -o data/L-X-A-kmers/L_profiling -p 1 -k 27 -n L_profiling
```

---

### Making tons of FigTree trees

It's a quick and dirty soltion to visualize all the tree toplogies for sanity checks

```bash
phylopath=data/GRC_phylogenies
for treefile in $phylopath/15_genetrees_L/*treefile; do
  pngfile=$(basename $treefile .treefile).png
  java -Xms64m -Xmx512m -jar lib/figtree.jar -graphic PNG $treefile $phylopath/plot_$pngfile
done
```

### More kmer explorations

Talbe of mapped kmers:

```
data/L-X-A-kmers/mapping/table_of_mapped_kmers_PacBio.tsv         # PacBio assembly
```

I will also need a list of lengths of PacBio contigs

```
samtools view -H data/L-X-A-kmers/mapping/A-27mer_mapped_racon6pe.bam | grep "^@SQ" | awk '{ print substr($2,4) "\t" substr($3,4) }' > data/L-X-A-kmers/mapping/table_of_mapped_kmers_spades.tsv
```

The quality of assignment for each is explored in [this script](scripts/kmer-assigment-of-L-X-A/exploring_mapped_kmers_in_assemblies.R). So it seems that the assignment is not perfect. Nearly all scaffolds have more than one category of kmers mapping on them. Usually it's not that bad, but in some cases it is. So manually inspect some of Illumina and PacBio assmeblies

```
echo -n -e {L,X,A}-27mer_mapped_spades_assembly.bam\\n | tr -d ' ' > illumina_bam_files.list
# illumina_bam_files.list do magic to get illumina_bam_files.bed
samtools depth -f illumina_bam_files.list -b illumina_scaffolds_to_inspect.bed > illumina_inspected_scfs.depth
samtools depth -f illumina_bam_files.list -r 'NODE_99_length_112107_cov_44.232839' > illumina_NODE_99_length_112107.depth
```

For PacBio

```
echo -n -e {L,X,A}-27mer_mapped_racon6pe.bam\\n | tr -d ' ' > racon6pe_bam_files.list
samtools depth -f racon6pe_bam_files.list -b  racon6pe_scaffolds_to_inspect.bed > racon6pe_inspected_scfs.depth
```

Check coverage depth of `ctg11` in PB asm it has practically no kmers mapping.

```
python3 scripts/kmer-assigment-of-L-X-A/kmer_depth2blockwise_depth.py data/L-X-A-kmers/mapping/illumina_inspected_scfs.depth data/L-X-A-kmers/mapping/table_of_mapped_kmers_spades.tsv
python3 scripts/kmer-assigment-of-L-X-A/kmer_depth2blockwise_depth.py data/L-X-A-kmers/mapping/illumina_NODE_99_length_112107.depth 500 > data/L-X-A-kmers/mapping/table_of_mapped_kmers_NODE_99_length_112107.tsv
# data/L-X-A-kmers/mapping/racon6pe_inspected_scfs.depth
```

and explore the tables in [scripts/kmer-assigment-of-L-X-A/plot_blockwise_kmer_assignment.R](scripts/kmer-assigment-of-L-X-A/plot_blockwise_kmer_assignment.R).

Ha, `illumina_NODE_99_length_112107` is a super clear chimera of X and L.

![NODE_99_chimera](https://user-images.githubusercontent.com/8181573/75560500-2e07f600-5a3d-11ea-81eb-baad5fd94646.png)