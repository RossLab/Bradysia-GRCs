### Sorting of Sciara genome

The goal is to figure reliable assignment of genome assembly of Scaira to X-linked (X), L-linked (L) and autosomal (A) sequences.
- [Project board](https://github.com/orgs/RossLab/projects/1).


Input data:
 - Illumina testes
 - Illumina head
 - Illumina testes assembly (spades)
 - PacBio testes
 - PacBio heads (with the X inversion)
 - PacBio testes polished assembly (ReadBeans, Racon with both long and short reads)

Using 2d kmer spectra of the Illumina testes and head samples we can quite well identify individual clouds of kmers that belong to A/X/L respectively. The first step will be to isolate the kmers. Then there are several options:
1. Categorise scaffolds of the Illumina asm; then map Illumina asm to PacBio asm
2. Categorise directly ctigs of the PacBio asm
3. Categorise PacBio reads and assemble separatelly A/X/L
4. Categoriaw Illumina reads and map the reads to improve the mapping (and the exact matches will make finishing out the reads much easier)

We aim for clear cut separation. If there will be any discrepancy in the signal, we need to understand why.

#### Isolation of kmers

Isolation of kmers that correspond to X, L or autosomes respectively. For 27-mers the coverage thresholds could be

A
 - 125 < head   < 175
 - 80  < testes < 140

X
 - 50  < head   < 100
 - 60  < testes < 100

L
 - head < 5
 - 10 < testes # this might be too stringent


Practically it will require using KMC to dump the kmers in alphabetical order (specify L=5 as we don't need all the error kmers). Then to write a python parser that reads the two files simuntaneously and merges the coverage information as

```
kmer	testes	head
ATTCAT	54	72
AGCGGC  24	1
...
```

with this then it will be quite streightforward to subselect sets described above.

##### actually doing it

Extract dumps of kmers using [kmc](https://github.com/refresh-bio/KMC).

```
conda activate default_genomics
# build kmer db
qsub -cwd -N kmc_heads -V -pe smp64 20 -b yes 'L=5; U=175; SCRATCH=/scratch/$USER/$JOB_ID/; mkdir -p $SCRATCH/tmp data/L-X-A-kmers/kmer_db; kmc -k27 -t20 -m64 -ci$L -cs$U data/L-X-A-kmers/raw_reads/bamfilter.head2.clc.mar29.scop.head2.vs.scop.clc.sorted.bam.InIn.fq.gz $SCRATCH/head_kmer_counts $SCRATCH/tmp && mv $SCRATCH/head_kmer_counts* data/L-X-A-kmers/kmer_db'
# generate alphabetically sorted dump of kmers and their coverages
qsub -cwd -N kmc_dump_heads -V -pe smp64 20 -b yes 'L=5; U=174; SCRATCH=/scratch/$USER/$JOB_ID; mkdir -p $SCRATCH; kmc_tools transform data/L-X-A-kmers/kmer_db/head_kmer_counts -ci$L -cx$U dump -s $SCRATCH/head_k27.dump && mv $SCRATCH/head_k27.dump data/L-X-A-kmers/kmer_db'
```

and the same for testes

```
qsub -cwd -N kmc_testes -V -pe smp64 20 -b yes 'L=5; U=175; SCRATCH=/scratch/$USER/$JOB_ID/; mkdir -p $SCRATCH/tmp data/L-X-A-kmers/kmer_db; kmc -k27 -t20 -m64 -ci$L -cs$U data/L-X-A-kmers/raw_reads/bamfilter.testes.clc.mar29.scop.testes.vs.scop.clc.sorted.bam.InIn.fq.gz $SCRATCH/testes_kmer_counts $SCRATCH/tmp && mv $SCRATCH/testes_kmer_counts* data/L-X-A-kmers/kmer_db'
qsub -cwd -N kmc_dump_testes -V -pe smp64 20 -b yes 'L=5; U=174; SCRATCH=/scratch/$USER/$JOB_ID/; mkdir -p $SCRATCH; kmc_tools transform data/L-X-A-kmers/kmer_db/testes_kmer_counts -ci$L -cx$U dump -s $SCRATCH/testes_k27.dump && mv $SCRATCH/testes_k27.dump data/L-X-A-kmers/kmer_db'
```

merge the two kmer dumps into one. Using a python script.

```
qsub -cwd -N dump_merging -V -pe smp64 1 -b yes 'SCRATCH=/scratch/$USER/$JOB_ID/; mkdir -p $SCRATCH; ./scripts/kmer-assigment-of-L-X-A/merge_two_dumps.py data/L-X-A-kmers/kmer_db/head_k27.dump data/L-X-A-kmers/kmer_db/testes_k27.dump > $SCRATCH/merged_k27.dump && mv $SCRATCH/merged_k27.dump data/L-X-A-kmers/kmer_db/; rmdir $SCRATCH'
```

This will take approximatelly an hour and half. If that works out I can use simply the threshold up there to create fasta files with A/X/L kmers. Using the treshold above we generate kmer fasta files using another python sript

```
qsub -cwd -N sort_out_kmers -V -pe smp64 1 -b yes 'scripts/kmer-assigment-of-L-X-A/dump2fasta.py data/L-X-A-kmers/kmer_db/merged_k27.dump'
```

Just to be sure that I gut more less that same 2d histogram as Christina I'll just plot it in R, although I know that it is not the most efficient way. To reduce the time, I just extract the coverages

```
SCRATCH=/scratch/$USER/$JOB_ID/
awk '{print $2 "\t" $3}' data/L-X-A-kmers/kmer_db/merged_k27.dump > $SCRATCH/merged_k27_cov_only.dump
Rscript scripts/kmer-assigment-of-L-X-A/plot_2d_histogram.R $SCRATCH/merged_k27_cov_only.dump
rm $SCRATCH/merged_k27_cov_only.dump && rmdir $SCRATCH'
```

Ha, one more thing. Is it possible to get the L size using kmers?

```
scripts/kmer-assigment-of-L-X-A/dump2hist.py data/L-X-A-kmers/kmer_db/merged_k27.dump data/L-X-A-kmers/L_kmers_k27.hist
```

generates the histogram that I will feed to genomescope

```
genomescope.R -i data/L-X-A-kmers/L_kmers_k27.hist -o data/L-X-A-kmers/L_profiling -p 1 -k 27 -n L_profiling
```

#### Options 1 and 2 - matching the kmers to assemblies

Christina mapped the kmers to both assemblies (TODO add the commands) and using `scripts/kmer-assigment-of-L-X-A/bams2kmer_tab.py` script she got the table of number of L/X/A kmers assigned to each scaffold

```
data/L-X-A-kmers/mapping/table_of_mapped_kmers_spades.tsv  # Illumina assembly
data/L-X-A-kmers/mapping/table_of_mapped_kmers_PacBio.tsv         # PacBio assembly
```

I will also need a list of lengths of PacBio contifs

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

##### L assigments

We was that the L chromosome are the easiest to assign - the kmer mapping is really sweet, most of the messed up assignments were between X and A. We have lots and lots of obviously L scaffolds that are nearly fully covered by L kmers and no others.

We want to be conservative therefore we use a combined information of read mapping (coverage analysis) and kmers. This will allow to define a core set Ls (hopefully nearly all of them) and probably a smaller set that have a kmer or coverage support.

The input for the assignment is

```
data/table.covdiff.germ.soma.txt # coverage table
data/L-X-A-kmers/mapping/table_of_mapped_kmers_spades.tsv # mapped kmers table
```

and it is done in [scripts/kmer-assigment-of-L-X-A/L-assignment.R](scripts/kmer-assigment-of-L-X-A/L-assignment.R) script. Generates `data/scaffold_assignment_tab_full.tsv` with all assignments and all the decision making metadata.
