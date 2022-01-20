## Objective
- to determine what the coverage of the L chromosome genes are in the Illumina assembly.
- we want to know if some genes are at half the coverage of others and if so if those are the ones that have paralogs (i.e. if the paralogs have been assembled as two separate things).
- this will give us an idea of how similar the L chromosomes are to each other.


## Approach:
- have genome annotation file that calls gene regions in the Illumina assembly.
- need a file that has testes (germ) reads mapped to genome .
- filter out only the L chromosome genes and look at cov histogram

#### mapping the testes library to the masked illumina assembly used for the annotation.
genome: 19_hisat2/scop.spades2.min1kb.trimmed.fasta.masked
reads: bamfilter.testes.clc.mar29.scop.testes.vs.scop.clc.sorted.bam.InIn.fq.gz

bwa index 19_hisat2/scop.spades2.min1kb.trimmed.fasta.masked
command
```
qsub -o logs -e logs -cwd -N bwamem -V -pe smp64 24 -b yes 'bwa mem -t 24 -p 19_hisat2/scop.spades2.min1kb.trimmed.fasta.masked bamfilter.testes.clc.mar29.scop.testes.vs.scop.clc.sorted.bam.InIn.fq.gz | samtools view -b - > /scratch/chodson/germ.vs.spades.masked.bam && rsync -avhP --remove-source-files /scratch/chodson/germ.vs.spades.masked.bam ../Sciara-L-chromosome/data/germ.vs.spades.masked.bam && samtools flagstat ../Sciara-L-chromosome/data/germ.vs.spades.masked.bam > ../Sciara-L-chromosome/data/germ.vs.spades.masked.stats.txt'
```

#### trying bedtools to get cov of annotated genes:

the way bedtools coverage works seems to be that you put the file that you want the cov for as -a and the one to take the info from as -b. So in this case the gff3 file should be -a. I've also filtered the gff3 file for just the genes so I only get what I want cov values for (otherwise I'm not sure exactly how big the file would be

in Sciara-L-chromosome directory
```
grep "gene" data/genome/annotation.gff3 > data/genome/sciara.gene.ILannotatin.gff3
```
```
qsub -o logs -e logs -cwd -N sort -V -pe smp64 4 -b yes 'bedtools sort -i data/genome/sciara.gene.ILannotatin.gff3 > data/genome/sciara.gene.ILannotatin.sort.gff3'

qsub -o logs -e logs -cwd -N sort -V -pe smp64 8 -b yes 'samtools sort -@ 8 -o data/germ.vs.spades.maskedsort.bam -O BAM -T /scratch/chodson/ data/germ.vs.spades.masked.bam'

qsub -o logs -e logs -cwd -N sort -V -pe smp64 8 -b yes 'bedtools coverage -mean -a data/genome/sciara.gene.ILannotatin.sort.gff3 -b data/germ.vs.spades.maskedsort.bam > data/gene.cov.braker.annotation.tsv'
```

ok, done, it gave me back the gff file info and the mean cov across each gene (I think it's computed as the coverage at each position averaged over the whole gene). Now going to filter the L contigs out of the tsv file.

#### extracting just the L chromosome genes and making a histogram of their coverage
In R.
scripts/cov_L_chromosome_genes.R

get a histogram of coverages of each L gene:
[need to ask Kamil how to add image of histogram]

we find that the two peaks in the histogram are at 24.64183 and 30.29147 cov, which means that there isn't any evidence that some of the L genes on homologs have been assembled as one thing and some have been assembled as two things (we would then expect a peak at half the cov). However, it's possible that what we're seeing is that most of the L genes on homologous chromosomes are not assembled as the same thing (i.e. we have two completely separate L chromosomes with distant homology to each other). If this is the case we would expect one homolog in one of the coverage peaks and one in the other. We're trying to figure out if that's the case...

### Getting per-scf coverage

I found we don't have a clear coverage per scaffold. I guess we have something like that in the assignmnet table (there is number of reads and length of scaffolds for germ and soma libraries). I will try to get a new file using the same mapped file that was used for Figure 4 (showing that the paralogs are likely on two different Ls). Then i will compare them, hopefully they will be moreless the same.

```
qsub -o logs -e logs -cwd -N genome_cov -V -pe smp64 1 -b yes '~/generic_genomics/depth2depth_per_contig.py <(samtools depth data/germ.vs.spades.maskedsort.bam) > data/genome/depth_per_scf.tsv'
```

I will add the coverage one column to the scaffold assignmnet table (`data/scaffold_assignment_tab_full.tsv`)

```{R}
cov_tab <- read.table('data/genome/depth_per_scf.tsv', col.names = c('scf', 'depth'), stringsAsFactors = F)
scf_asn <- read.delim("data/scaffold_assignment_tab_full.tsv", header=T, stringsAsFactors = F)

row.names(cov_tab) <- cov_tab$scf
cov_tab$len <- NA
cov_tab[scf_asn$scf, 'len'] <- scf_asn$len
cov_tab$cov <- cov_tab$depth / cov_tab$len
row.names(scf_asn) <- scf_asn$scf

scf_asn$cov <- NA
scf_asn[scf_asn$scf, 'cov'] <- cov_tab$cov

write.table(scf_asn, 'data/scaffold_assignment_tab_full.tsv', quote = F, sep = "\t", row.names = F)
```

and now I can try to get putative L1 and L2 assignments (for clarity in an indipendent R session)

```R
scf_asn <- read.delim("data/scaffold_assignment_tab_full.tsv", header=T, stringsAsFactors = F)

L_scfs <- scf_asn[scf_asn$assignment == 'L', ]

hist(log10(L_scfs$cov), breaks = 300) # , xlim = c(1, 2)

# the two peaks around 1.4 and 1.5 are the L1 and L2 we are looking for. The other peaks represent paralogs or shared sequences between L1 and L2. Unfortunatelly, impossible to disentagle at this point

# the only assignable scaffolds are the single copy scaffolds, so let's subset them:
L_scfs <- L_scfs[L_scfs$cov < 40, ]

hist(L_scfs$cov, breaks = 100)
# this signal is much clearer now. It's very likely that loads of the noise come from short scaffolds

plot(L_scfs$cov, L_scfs$len)
lines(c(0, 100), c(15000, 15000), lty = 2, col = 'red')
# yeah, and their assignments would not be reliable anyway, so let's just assign all the scaffolds > 15,000

# why there are NA values for scf length? TODO: look into it!
L_scfs <- L_scfs[!is.na(L_scfs$len), ]

sum(L_scfs$len[L_scfs$len > 15000])
# Also, >15k is still 104M of the L scaffolds, that's solid for assignments

hist(L_scfs$cov, breaks = 100)

# finally, I will remove the few scaffolds with suspiciously low coverage (< 20)

L_scfs <- L_scfs[L_scfs$len > 15000 & L_scfs$cov > 15, ]

sum(L_scfs$len) # still 104M (it was only a few kb removed)
hist(L_scfs$cov, breaks = 100)
```

#### mixed models

```R
library("mixtools")
library("ggplot2")
library("magrittr")
library('plotGMM')

scf_asn <- read.delim("data/scaffold_assignment_tab_full.tsv", header=T, stringsAsFactors = F)

L_scfs <- scf_asn[scf_asn$assignment == 'L', ]
L_scfs <- L_scfs[!is.na(L_scfs$len), ]
L_scfs <- L_scfs[L_scfs$len > 15000 & L_scfs$cov > 17 & L_scfs$cov < 37, ]

mixmdl_k2 <- normalmixEM(L_scfs$cov, k = 2, mu = c(25, 31))

data.frame(x = mixmdl_k2$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = 0.5, colour = "black", alpha=1) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl_k2$mu[1], mixmdl_k2$sigma[1], lam = mixmdl_k2$lambda[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl_k2$mu[2], mixmdl_k2$sigma[2], lam = mixmdl_k2$lambda[2]),
                colour = "blue", lwd = 1.5) +
  geom_vline(aes(xintercept=mixmdl_k2$mu[1]),linetype="dashed", size=1.5,colour="red")+
  geom_vline(aes(xintercept=mixmdl_k2$mu[2]),linetype="dashed",size=1.5,colour="blue")+
  theme_bw()+
  ggtitle("Allacma fusca")+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        plot.title = element_text(size=20))+
  xlab("Coverage")+
  ylab("Density")

# TODO using mu and sd, calculate propabilities of individual scaffolds being from that distribution

# L1 (1 - ) because it's only the right tail of the distribution we consider
L1_asn_prob <- 1 - pnorm(L_scfs$cov, mixmdl_k2$mu[1], mixmdl_k2$sigma[1])
# L2
L2_asn_prob <- pnorm(L_scfs$cov, mixmdl_k2$mu[2], mixmdl_k2$sigma[2])

L_scfs$L_asn <- "L2"
L_scfs$L_asn[L1_asn_prob > L2_asn_prob] <- "L1"

L1s <- L_scfs$L_asn == "L1"

L_scfs$L_asn_LL <- NA
# adding "+ 1e6" to denominator so the log values are never NA
L_scfs[L1s, 'L_asn_LL'] <- log(L1_asn_prob[L1s] / (L2_asn_prob[L1s] + 1e6))
L_scfs[!L1s, 'L_asn_LL'] <- log(L2_asn_prob[!L1s] / (L1_asn_prob[!L1s] + 1e6))

plot(L_scfs$L_asn_LL ~ L_scfs$cov)
# Finally, it makes sense. Dropout of likelihoods around the divide of the two distributions

hist(L_scfs$L_asn_LL, breaks = 120)
# finally, removing values with likelihood smaller than 1e-16 -> those were too close to the divide of the two distribution for reliable assignment
#
L_scfs[L_scfs$L_asn_LL < -16, "L_asn"] <- paste0(L_scfs[L_scfs$L_asn_LL < -16, "L_asn"], 'c')

sum(L_scfs[L_scfs$L_asn == "L1", 'len']) / 1e6
# 34M of L1
# vs
sum(L_scfs[L_scfs$L_asn == "L2", 'len']) / 1e6
# 56M of L2
# probably means L2 will have some false positives in there as we would expect approximatelly the same size of single copy lengths of one and the other L chromosome

out_table <- L_scfs[, c('scf', 'len', 'cov', 'L_asn', 'L_asn_LL')]

write.table(out_table, 'data/L_scaffolds_chrom_assignment.tsv', quote = F, sep = "\t", row.names = F)
```

So, file `data/L_scaffolds_chrom_assignment.tsv` contains 2459 scaffolds wirh their assignmnets and corresponding likelihoods.

---

TODO: Recheck with Christina how mapping was done for 1 - assinging X/A/L and 2. getting coverage. They are not 100% correlated

See

```R
scf_asn <- read.delim("data/scaffold_assignment_tab_full.tsv", header=T, stringsAsFactors = F)

plot(150 * (testes.readcount / len) ~ cov, data = scf_asn, xlim = c(0, 300), ylim = c(0, 300))
```
