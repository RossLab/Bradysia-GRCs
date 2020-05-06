## goal: to determine what the coverage of the L chromosome genes are in the Illumina assembly. 
- we want to know if some genes are at half the coverage of others and if so if those are the ones that have paralogs (i.e. if the paralogs have been assembled as two separate things).
- this will give us an idea of how similar the L chromosomes are to each other. 


## Approach:
- have genome annotation file that calls gene regions in the Illumina assembly.
- need a file that has testes (germ) reads mapped to genome (already have for unmasked assembly but I think I'll redo with the masked assembly so I know I have the right thing).
- need to find a tool that uses this info to calculate the coverage for each gene. Can use samtools mpileup (but need to adjust so it only gives you coverage for genes not every base) but I'm going to try bedtools coverage since that says it will give cov with just bam file and gff annotation file (we'll see how this goes).
- filter out only the L chromosome genes and look at cov histogram

#### mapping the testes library to the masked illumina assembly used for the annotation.
genome: 19_hisat2/scop.spades2.min1kb.trimmed.fasta.masked
reads: bamfilter.testes.clc.mar29.scop.testes.vs.scop.clc.sorted.bam.ExIn.fq.gz

bwa index 19_hisat2/scop.spades2.min1kb.trimmed.fasta.masked
command
```
qsub -o logs -e logs -cwd -N bwamem -V -pe smp64 24 -b yes 'bwa mem -t 24 -p 19_hisat2/scop.spades2.min1kb.trimmed.fasta.masked bamfilter.testes.clc.mar29.scop.testes.vs.scop.clc.sorted.bam.InIn.fq.gz | samtools view -b - > /scratch/chodson/germ.vs.spades.masked.bam && rsync -avhP --remove-source-files /scratch/chodson/germ.vs.spades.masked.bam ../Sciara-L-chromosome/data/germ.vs.spades.masked.bam && samtools flagstat ../Sciara-L-chromosome/data/germ.vs.spades.masked.bam > ../Sciara-L-chromosome/data/germ.vs.spades.masked.stats.txt'
```

download bedtools
on chodson
```
conda install -c bioconda bedtools
```

ex. command
bedtools coverage [OPTIONS] -a <FILE> \
                             -b <FILE1, FILE2, ..., FILEN>
                             
other options for this step:
samtools mpileup

with ...

def getReference(refFasta, chrom, start, end):
    region = "%s:%d-%d" % (chrom, start, end)
    samtoolsOut = check_output(["samtools", "faidx", refFasta, region])
    # sys.stderr.write("".join([samtools, "faidx" + refFasta + region]))
    # sys.stderr.write(str(samtoolsOut))
    refSeq = ""
    for seq in str(samtoolsOut).split('\n'):
        if not seq.startswith(">"):
            refSeq += seq

    return refSeq.upper()


#### trying bedtools to get cov of annotated genes:

the way bedtools coverage works seems to be that you put the file that you want the cov for as -a and the one to take the info from as -b. So in this case the gff3 file should be -a. I've also filtered the gff3 file for just the genes so I only get what I want cov values for (otherwise I'm not sure exactly how big the file would be

in Sciara-L-chromosome directory
```
grep "gene" data/genome/annotation.gff3 > data/genome/sciara.gene.ILannotatin.gff3
```
I think I might sort some of the files before I run the cov command
(also, the samtools in chodson is acting weird, I'm using rnaalign env right now)

-hist (computes histogram for each feature in A. Don't want hist for all genes so I'm going to see what I get without this first)
```
qsub -o logs -e logs -cwd -N sort -V -pe smp64 4 -b yes 'bedtools sort -i data/genome/sciara.gene.ILannotatin.gff3 > data/genome/sciara.gene.ILannotatin.sort.gff3'

qsub -o logs -e logs -cwd -N sort -V -pe smp64 8 -b yes 'samtools sort -@ 8 -o data/germ.vs.spades.maskedsort.bam -O BAM -T /scratch/chodson/ data/germ.vs.spades.masked.bam'

qsub -o logs -e logs -cwd -N sort -V -pe smp64 8 -b yes 'bedtools coverage -mean -a data/genome/sciara.gene.ILannotatin.sort.gff3 -b data/germ.vs.spades.maskedsort.bam > data/gene.cov.braker.annotation.tsv'
```

ok, done, it gave me back the gff file info and the mean cov across each gene (I think it's computed as the coverage at each position averaged over the whole gene. Now going to filter the L contigs out of the tsv file.

#### extracting just the L chromosome genes and making a histogram of their coverage
In R. 
scripts/cov_L_chromosome_genes.R

get a histogram of coverages of each L gene:
[need to ask Kamil how to add image of histogram]

we find that the two peaks in the histogram are at 24.64183 and 30.29147 cov, which means that there isn't any evidence that some of the L genes on homologs have been assembled as one thing and some have been assembled as two things (we would then expect a peak at half the cov). However, it's possible that what we're seeing is that most of the L genes on homologous chromosomes are not assembled as the same thing (i.e. we have two completely separate L chromosomes with distant paralogy to each other). If this is the case we would expect one paralog in one of the coverage peaks and one in the other. We're trying to figure out if that's the case...
