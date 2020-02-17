### The assembly of Sciara

This is an exprimental approach for an assembly of the fungus gnat genome using combined information of two PacBio and Illumina datasets. The progress is recorded [here](https://github.com/RossLab/projects/1).

#### Step 1 - isolation of kmers

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
qsub -cwd -N kmc_heads -V -pe smp64 20 -b yes 'L=5; U=175; SCRATCH=/scratch/$USER/$JOB_ID/; mkdir -p $SCRATCH/tmp data/kmer_db; kmc -k27 -t20 -m64 -ci$L -cs$U data/raw_reads/bamfilter.head2.clc.mar29.scop.head2.vs.scop.clc.sorted.bam.InIn.fq.gz $SCRATCH/head_kmer_counts $SCRATCH/tmp && mv $SCRATCH/head_kmer_counts* data/kmer_db'
# generate alphabetically sorted dump of kmers and their coverages
qsub -cwd -N kmc_dump_heads -V -pe smp64 20 -b yes 'L=5; U=174; SCRATCH=/scratch/$USER/$JOB_ID; mkdir -p $SCRATCH; kmc_tools transform data/kmer_db/head_kmer_counts -ci$L -cx$U dump -s $SCRATCH/head_k27.dump && mv $SCRATCH/head_k27.dump data/kmer_db'
```

and the same for testes

```
qsub -cwd -N kmc_testes -V -pe smp64 20 -b yes 'L=5; U=175; SCRATCH=/scratch/$USER/$JOB_ID/; mkdir -p $SCRATCH/tmp data/kmer_db; kmc -k27 -t20 -m64 -ci$L -cs$U data/raw_reads/bamfilter.testes.clc.mar29.scop.testes.vs.scop.clc.sorted.bam.InIn.fq.gz $SCRATCH/testes_kmer_counts $SCRATCH/tmp && mv $SCRATCH/testes_kmer_counts* data/kmer_db'
qsub -cwd -N kmc_dump_testes -V -pe smp64 20 -b yes 'L=5; U=174; SCRATCH=/scratch/$USER/$JOB_ID/; mkdir -p $SCRATCH; kmc_tools transform data/kmer_db/testes_kmer_counts -ci$L -cx$U dump -s $SCRATCH/testes_k27.dump && mv $SCRATCH/testes_k27.dump data/kmer_db'
```

merge the two kmer dumps into one. Using a python script.

```
qsub -cwd -N dump_merging -V -pe smp64 1 -b yes 'SCRATCH=/scratch/$USER/$JOB_ID/; mkdir -p $SCRATCH; ./scripts/merge_two_dumps.py data/kmer_db/head_k27.dump data/kmer_db/testes_k27.dump > $SCRATCH/merged_k27.dump && mv $SCRATCH/merged_k27.dump data/kmer_db/; rmdir $SCRATCH'
```

This will take approximatelly an hour and half. If that works out I can use simply the threshold up there to create fasta files with A/X/L kmers. Using the treshold above we generate kmer fasta files using another python sript

```
qsub -cwd -N sort_out_kmers -V -pe smp64 1 -b yes 'scripts/dump2fasta.py data/kmer_db/merged_k27.dump'
```

Just to be sure that I gut more less that same 2d histogram as Christina I'll just plot it in R, although I know that it is not the most efficient way. To reduce the time, I just extract the coverages

```
SCRATCH=/scratch/$USER/$JOB_ID/
awk '{print $2 "\t" $3}' data/kmer_db/merged_k27.dump > $SCRATCH/merged_k27_cov_only.dump 
Rscript scripts/plot_2d_histogram.R $SCRATCH/merged_k27_cov_only.dump 
rm $SCRATCH/merged_k27_cov_only.dump && rmdir $SCRATCH'
```

#### Step 2 - matching the kmers to long reads

KAT has a tool that subselect sequences that have a kmer in a hash. I could build hashes out of the kmer subsets and then use that function to subselect PB reads. However, this won't allow a exploring.

##### Exploring mappability of kmers

Before I will attempt a scalable solution. I will attempt to take a relatively small set of reads (5000?) and map the three kemrs on them. That should be relatively fast (actualy that should be nearly that same fast as mapping to all the reads, but the difference is that I won't have to deal with memory).

I will try a `bwa mem` quick and dirty approach suggested [here](https://bioinformatics.stackexchange.com/a/7299/57).

```
zcat /data/ross/sequencing/raw/cgr.liv.ac.uk/pbio/LIMS21670_c553b0c0ac9a6d5a/1/FilteredSubreads/subreads.fastq.gz | sed -n '1~4s/^@/>/p;2~4p' | head -10000 > data/raw_reads/PB_sample.fasta
bwa index data/raw_reads/PB_sample.fasta
# data/kmers_k27_X.fasta
# data/kmers_k27_A.fasta
# data/kmers_k27_L.fasta
qsub -cwd -N map_X_kmers -V -pe smp64 16 -b yes 'bwa mem -t 10 -k 27 -T 27 -a -c 5000 data/raw_reads/PB_sample.fasta data/kmers_k27_X.fasta | samtools sort -@6 -O bam - > data/X-27mer_mapped_to_sample.bam'
qsub -cwd -N map_L_kmers -V -pe smp64 16 -b yes 'bwa mem -t 10 -k 27 -T 27 -a -c 5000 data/raw_reads/PB_sample.fasta data/kmers_k27_L.fasta | samtools sort -@6 -O bam - > data/L-27mer_mapped_to_sample.bam'
qsub -cwd -N map_A_kmers -V -pe smp64 16 -b yes 'bwa mem -t 10 -k 27 -T 27 -a -c 5000 data/raw_reads/PB_sample.fasta data/kmers_k27_A.fasta | samtools sort -@6 -O bam - > data/A-27mer_mapped_to_sample.bam'
```


Let's try to build bwa index on the full read set.

```
zcat /data/ross/sequencing/raw/cgr.liv.ac.uk/pbio/LIMS21670_c553b0c0ac9a6d5a/1/FilteredSubreads/subreads.fastq.gz | sed -n '1~4s/^@/>/p;2~4p' > data/raw_reads/PB_testes.fasta
qsub -cwd -N build_bwa -V -pe smp64 1 -b yes 'bwa index data/raw_reads/PB_testes.fasta'
```

and map the smallest kmer set to them

```
qsub -cwd -N map_X_kmers -V -pe smp64 32 -b yes 'bwa mem -t 26 -k 27 -T 27 -a -c 5000 data/raw_reads/PB_testes.fasta data/kmers_k27_X.fasta | samtools sort -@6 -O bam - > data/mapping_testes_27mer_X.bam'
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
