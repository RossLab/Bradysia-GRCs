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

Finally, merge the two kmer dumps into one. Using a python script.

```
# local testing
./scripts/merge_two_dumps.py data/kmer_db/head_k27_sample.dump data/kmer_db/testes_k27_sample.dump > data/kmer_db/merged_k27_sample.dump
```

```{Rscript}
library(hexbin)
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

kmer_dump <- read.table('data/kmer_db/merged_k27_sample.dump', col.names = c('kmer', 'heads', 'testes'))
hexbinplot(testes ~ heads, data=kmer_dump, colramp=rf)
```

```
qsub -cwd -N dump_merging -V -pe smp64 1 -b yes 'SCRATCH=/scratch/$USER/$JOB_ID/; mkdir -p $SCRATCH; ./scripts/merge_two_dumps.py data/kmer_db/head_k27.dump data/kmer_db/testes_k27.dump > $SCRATCH/merged_k27.dump && mv $SCRATCH/merged_k27.dump data/kmer_db/; rmdir $SCRATCH'
```

#### Step 2 - matching the kmers to long reads

KAT has a tool that subselect sequences that have a kmer in a hash. I could build hashes out of the kmer subsets and then use that function to subselect PB reads. However, this won't allow a exploring.

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
