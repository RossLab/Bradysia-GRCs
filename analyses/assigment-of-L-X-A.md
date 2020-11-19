## Sorting of Sciara genome

The goal is to figure reliable assignment of genome assembly of Scaira to X-linked (X), L-linked (L) and autosomal (A) sequences.
- [Project board](https://github.com/orgs/RossLab/projects/1).


Input data:
 - Illumina testes
 - Illumina head
 - Illumina testes assembly (spades)
 - PacBio testes
 - PacBio heads (with the X inversion)
 - PacBio testes polished assembly (ReadBeans, Racon with both long and short reads)

We use two approaches, the coverage ratio of head and testes libraries and using 2d kmer spectra.

### Table of content

- [Coverage analysis](#coverage-analysis)
  * [1. make/get bam file in which reads are mapped to contigs in assembly](#1-make-get-bam-file-in-which-reads-are-mapped-to-contigs-in-assembly)
  * [2. count mapped reads for each contig (and sorting bam files)](#2-count-mapped-reads-for-each-contig--and-sorting-bam-files-)
  * [3./4. compare reads from testes to head library of log2 scale and make histogram](#3-4-compare-reads-from-testes-to-head-library-of-log2-scale-and-make-histogram)
  * [Generating some summary stats for the bam files.](#generating-some-summary-stats-for-the-bam-files)
- [Kmer spectra based assigment analyses](#kmer-spectra-based-assigment-analyses)
  * [Isolation of kmers](#isolation-of-kmers)
  * [actually doing it](#actually-doing-it)
  * [Matching the kmers to assemblies](#options-1-and-2---matching-the-kmers-to-assemblies)
- [Final L assigments](#final-l-assigments)

### Coverage analysis

The purpose of this script is to map the reads of different libraries (in this case germline and somatic) back to a reference genome and then look for places where there are differences in how many reads map back to the genome from each library to identify regions at different copy numbers in the different libraries. Here's a list of the steps:
1. map reads from each library to the reference genome (bwa mem)
2. count the number of reads mapping to each contig from each library ()
3. import the counts into R and compare the log2 difference between the libraries. (you do the log2 difference so the values are centered at 0 if the contig is at equal coverage in both libraries, 1 for things at double the coverage in the first library (for example, a contigs from a chromosome with two copies in the first tissue and 1 in the second tissue type), and -1 for things at half the coverage in the first library compared to the second).
4. then look at the histogram of the log2 difference and decide where to set cutoffs for different regions of the genome (it's good to have some sort of expectation before hand for where your cutoffs are likely to be but you should also see what the histogram actually looks like).
5. then apply the cutoffs to the data and assign different contigs to different regions of the genome (I generated a table with this information). The most important thing is to have the assignments and the contig ID's in the same table.

Files

```
/data/ross/mealybugs/analyses/sciara_coprophila/12_Lcontigs/scop.testeslib.vs.spades2.may8.sam
/data/ross/mealybugs/analyses/sciara_coprophila/12_Lcontigs/scop.headlib.vs.spades2.may8.sam
```

#### Here are some important things about what I did below
- used the min1kb assembly for mapping

#### 1. make/get bam file in which reads are mapped to contigs in assembly

```
#command to change from sam -S part of command to bam -b part of command
samtools view -S -b sample.sam > sample.bam
# to view file...
samtools view sample.bam | head
```

mapping command:
```
bwa index 10_spades2/contigs.min1kb.fasta && bwa mem -t 16 -p 10_spades2/contigs.min1kb.fasta bamfilter.head2.clc.mar29.scop.head2.vs.scop.clc.sorted.bam.InIn.fq.gz | samtools view -b - > /scratch/chodson/scop.head2.vs.spades2.bam && rsync --remove-source-files /scratch/chodson/scop.head2.vs.spades2.bam 10_spades2/ && touch 10_spades2/bwa.head.done

bwa mem -t 16 -p 10_spades2/contigs.min1kb.fasta bamfilter.testes.clc.mar29.scop.testes.vs.scop.clc.sorted.bam.InIn.fq.gz | samtools view -b - > /scratch/chodson/scop.testes.vs.spades2.bam && rsync --remove-source-files /scratch/chodson/scop.testes.vs.spades2.bam 10_spades2/ && touch 10_spades2/bwa.testes.done
```

#### 2. count mapped reads for each contig (and sorting bam files)

e.x. command

```
samtools sort -o male_sort.bam male_to_mf.bam ## sort bam
samtools index male_sort.bam ## index
samtools idxstats male_sort.bam > male_reads_mapped_unique.txt ## output mapped read counts
```

my command:
```
samtools sort -@ 8 -o 13_covanalysis/head_sort.bam 10_spades2/scop.head2.vs.spades2.bam && samtools index 13_covanalysis/head_sort.bam && samtools idxstats 13_covanalysis/head_sort.bam > 13_covanalysis/head_reads_mapped_unique.txt
samtools sort -@ 8 -o 13_covanalysis/testes_sort.bam 10_spades2/scop.testes.vs.spades2.bam && samtools index 13_covanalysis/testes_sort.bam && samtools idxstats 13_covanalysis/testes_sort.bam > 13_covanalysis/testes_reads_mapped_unique.txt
```

#### 3./4. compare reads from testes to head library of log2 scale and make histogram

- all this was done by moving count files to local computer and running in R

```R
# reading in files
setwd("/Users//christina//Dropbox//Sciara//assembly.outputs//covanalysis")
testes.read.counts<- read.delim("testes_reads_mapped_unique.txt",header=FALSE)
head.read.counts<- read.delim("head_reads_mapped_unique.txt",header=FALSE)

# log(x/y) = log(x) - log(y)

# getting cov dif measurement between libraries
logdif2<-log2(testes.read.counts[ ,3])-log2(head.read.counts[ ,3])
cbind(testes.read.counts,logdif)

# trying out some histograms
hist(logdif2,breaks=500)
hist(logdif2,breaks=200, ylim=c(0, 3000),xlim=c(-1, 8))
hist(logdif2,breaks=500, xlim=c(-3, 3))

# making an output file with all columns of interest
ID<-testes.read.counts[ ,1]
seqlength<-testes.read.counts[ ,2]
testes.readcount<-testes.read.counts[ ,3]
head.readcount<-head.read.counts[ ,3]
testes.head.cov.diff<-data.frame(ID,seqlength,testes.readcount, head.readcount,logdif,logdif2)
head(testes.head.cov.diff)


write.table(testes.head.cov.diff, file='testes.head.cov.diff.tsv', quote=FALSE, sep='\t', col.names = NA)
###looks like above worked pretty well
```

Generating file with just the contigs at higher cov in germ tissue (log difference of greater than 2) and with a length greater than 1000

```
cat testes.head.cov.diff.tsv | awk '($7> 2)' > highcov.germ.txt
cat highcov.germ.txt | awk '($3> 1000)'> highcov.germ2.txt

wc -l highcov.germ.txt highcov.germ2.txt
 8441 highcov.germ.txt
 8438 highcov.germ2.txt
```

#### Generating some summary stats for the bam files.
```
samtools flagstat 13_covanalysis/head_sort.bam > 13_covanalysis/head_bam_stats.txt && samtools flagstat 13_covanalysis/testes_sort.bam > 13_covanalysis/testes_bam_stats.txt
```

For `head_sort.bam`:

```
341502703 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
12275879 + 0 supplementary
0 + 0 duplicates
334865215 + 0 mapped (98.06% : N/A)
329226824 + 0 paired in sequencing
164613412 + 0 read1
164613412 + 0 read2
261651594 + 0 properly paired (79.47% : N/A)
318858184 + 0 with itself and mate mapped
3731152 + 0 singletons (1.13% : N/A)
56526900 + 0 with mate mapped to a different chr
40249475 + 0 with mate mapped to a different chr (mapQ>=5)
```

For `testes_sort.bam`:

```
324311638 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
12501284 + 0 supplementary
0 + 0 duplicates
318275962 + 0 mapped (98.14% : N/A)
311810354 + 0 paired in sequencing
155905177 + 0 read1
155905177 + 0 read2
248931976 + 0 properly paired (79.83% : N/A)
302354384 + 0 with itself and mate mapped
3420294 + 0 singletons (1.10% : N/A)
52747120 + 0 with mate mapped to a different chr
38122814 + 0 with mate mapped to a different chr (mapQ>=5)
```

### Kmer spectra based assigment analyses

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

#### actually doing it

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

#### Matching the kmers to assemblies

Christina mapped the kmers to both assemblies (`bwa mem -k 27 -T 27 -a -c 5000`) and using `scripts/kmer-assigment-of-L-X-A/bams2kmer_tab.py` script she got the table of number of L/X/A kmers assigned to each scaffold

```
data/L-X-A-kmers/mapping/table_of_mapped_kmers_spades.tsv  # Illumina assembly
```

The stats of mapped kmers are shown on following plots (generated by `scripts/kmer-assigment-of-L-X-A/L-assignment_plots.R` script):

TODO

TODO

TODO

TODO

### Final L assigments

We was that the L chromosome are the easiest to assign - the kmer mapping is really sweet, most of the messed up assignments were between X and A. We have lots and lots of obviously L scaffolds that are nearly fully covered by L kmers and no others.

We want to be conservative therefore we use a combined information of read mapping (coverage analysis) and kmers. This will allow to define a core set Ls (hopefully nearly all of them) and probably a smaller set that have a kmer or coverage support.

The input for the assignment is

```
data/table.covdiff.germ.soma.txt # coverage table
data/L-X-A-kmers/mapping/table_of_mapped_kmers_spades.tsv # mapped kmers table
```

and it is done in [scripts/kmer-assigment-of-L-X-A/L-assignment.R](scripts/kmer-assigment-of-L-X-A/L-assignment.R) script. Generates `data/scaffold_assignment_tab_full.tsv` with all assignments and all the decision making metadata.
