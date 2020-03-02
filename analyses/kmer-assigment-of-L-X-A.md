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
data/highcov_germ.tsv # coverage table
data/L-X-A-kmers/mapping/table_of_mapped_kmers_spades.tsv # mapped kmers table
```

and it is done in [scripts/kmer-assigment-of-L-X-A/L-assignment.R](scripts/kmer-assigment-of-L-X-A/L-assignment.R) script.

#### Option 3 - matching the kmers to long reads

KAT has a tool that subselect sequences that have a kmer in a hash. I could build hashes out of the kmer subsets and then use that function to subselect PB reads. However, this won't allow a exploring.

##### Exploring mappability of kmers

Before I will attempt a scalable solution. I will attempt to take a relatively small set of reads (5000?) and map the three kemrs on them. That should be relatively fast (actualy that should be nearly that same fast as mapping to all the reads, but the difference is that I won't have to deal with memory).

I will try a `bwa mem` quick and dirty approach suggested [here](https://bioinformatics.stackexchange.com/a/7299/57).

```
zcat /data/L-X-A-kmers/ross/sequencing/raw/cgr.liv.ac.uk/pbio/LIMS21670_c553b0c0ac9a6d5a/1/FilteredSubreads/subreads.fastq.gz | sed -n '1~4s/^@/>/p;2~4p' | head -10000 > data/L-X-A-kmers/raw_reads/PB_sample.fasta
bwa index data/L-X-A-kmers/raw_reads/PB_sample.fasta
# data/L-X-A-kmers/kmers_k27_X.fasta
# data/L-X-A-kmers/kmers_k27_A.fasta
# data/L-X-A-kmers/kmers_k27_L.fasta
qsub -cwd -N map_X_kmers -V -pe smp64 16 -b yes 'bwa mem -t 10 -k 27 -T 27 -a -c 5000 data/L-X-A-kmers/raw_reads/PB_sample.fasta data/L-X-A-kmers/kmers_k27_X.fasta | samtools sort -@6 -O bam - > data/L-X-A-kmers/X-27mer_mapped_to_sample.bam'
qsub -cwd -N map_L_kmers -V -pe smp64 16 -b yes 'bwa mem -t 10 -k 27 -T 27 -a -c 5000 data/L-X-A-kmers/raw_reads/PB_sample.fasta data/L-X-A-kmers/kmers_k27_L.fasta | samtools sort -@6 -O bam - > data/L-X-A-kmers/L-27mer_mapped_to_sample.bam'
qsub -cwd -N map_A_kmers -V -pe smp64 16 -b yes 'bwa mem -t 10 -k 27 -T 27 -a -c 5000 data/L-X-A-kmers/raw_reads/PB_sample.fasta data/L-X-A-kmers/kmers_k27_A.fasta | samtools sort -@6 -O bam - > data/L-X-A-kmers/A-27mer_mapped_to_sample.bam'
```


Let's try to build bwa index on the full read set.

```
zcat /data/L-X-A-kmers/ross/sequencing/raw/cgr.liv.ac.uk/pbio/LIMS21670_c553b0c0ac9a6d5a/1/FilteredSubreads/subreads.fastq.gz | sed -n '1~4s/^@/>/p;2~4p' > data/L-X-A-kmers/raw_reads/PB_testes.fasta
qsub -cwd -N build_bwa -V -pe smp64 1 -b yes 'bwa index data/L-X-A-kmers/raw_reads/PB_testes.fasta'
```

and map the smallest kmer set to them

```
qsub -cwd -N map_X_kmers -V -pe smp64 32 -b yes 'bwa mem -t 26 -k 27 -T 27 -a -c 5000 data/L-X-A-kmers/raw_reads/PB_testes.fasta data/L-X-A-kmers/kmers_k27_X.fasta | samtools sort -@6 -O bam - > data/L-X-A-kmers/mapping_testes_27mer_X.bam'
```

### Potential tweaks

I do see some potential pitfalls, here I write what we could do to solve them

#### lower the k

If we will have troubles to find enough matches of kmers in long reads we could lower the k. There will be less kmers to work with, but the probability of a kmer being correctly called in the PacBio read will increase dramatically.

Here for different kmer lengths probabilities of a read kmer matching representing the correct genomic sequence assuming 9% error rate in sequencing.

k = 27; 0.078;
k = 21; 0.138;
k = 17; 0.201;

7.8% of kmers matching seems fine to me even for the shorter reads in the dataset. For sanity reasons I will do some simple simulations

```{R}
# function to run one replicate
simulate_read <- function(k, error_rate, read_length){ simulated_bases <- rbinom(n=read_length, prob=1-error_rate, size = 1); kmers = 0; for(i in 1:(read_length - k)){if(all(simulated_bases[i:(i+k)])){ kmers = kmers + 1 }}; kmers }
# running 1000 replicates
k27_rl1000 <- sapply(rep(27, 1000), simulate_read, 0.09, 1000)
quantile(k27_rl1000, 0.01)
> 1%
> 11
```

i.e. 99% of reads at least 1000 bases long should have more than 11 kmers matching, which should be more than sufficient to classify the read. However, here I don't take into account that a certain proportion of the genome is not in a single copy.
