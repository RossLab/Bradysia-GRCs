

#setwd("/Users/christina/projects/Sciara-L-chromosome/")
library(tximport)
scf_asn<- read.delim("data/scaffold_assignment_tab_full.tsv",header=T, stringsAsFactors = F)
trans_asn<- read.delim("data/salmon_rnaseq/transcript.scaffold.map.tsv",header=FALSE, stringsAsFactors = F)
scf_asn[is.na(scf_asn)] <- "NA"

trans_asn$transcript<- trans_asn$V1
trans_asn$scf<- trans_asn$V2
full_asn <- merge(scf_asn, trans_asn)

#table that has transcript ID's and scaffold name and assignment
transcript_asn <- full_asn[ c(1,2,3,6,11,15,18) ]

#small summary fig of number of transcripts belonging to each chromosome type
transcripts.per.chromosome <- table(transcript_asn$assignments)[c(1:4,6:7,5)]
pal <- c('green','green', 'yellow','yellow','red', 'red', 'grey')
barplot(transcripts.per.chromosome, col=pal, xlab='chromosome assignment', ylab='number of transcripts')

##reading in rnaseq data(the paths to this data might need to be changed on qm (it's currently in sciara_coprophila/20_braker/ I think))
samples <- read.table('data/salmon_rnaseq/RNA_sample_info.tsv', header=T, stringsAsFactors = F)
samples$filenames <- file.path("data/salmon_rnaseq/1_expression/", samples$sample, "quant.sf")


tx2gene <- read.table('data/salmon_rnaseq/scop.braker_annot1.IL.transcripts2genes.map', col.names = c('GENEID', 'TXNAME'), stringsAsFactors = F)[,c(2,1)]
txi <- tximport(samples$filenames, type = "salmon", tx2gene = tx2gene)

row.names(tx2gene) <- tx2gene$TXNAME
L_transcripts <- tx2gene[transcript_asn[transcript_asn$assignments == 'L', 'transcript'], 'GENEID']
X_transcripts <- tx2gene[transcript_asn[transcript_asn$assignments == 'X', 'transcript'], 'GENEID']
A_transcripts <- tx2gene[transcript_asn[transcript_asn$assignments == 'A', 'transcript'], 'GENEID']

exp_matrix <- txi[[1]]

L_transcripts_quantified <- L_transcripts[L_transcripts %in% row.names(exp_matrix)]
X_transcripts_quantified <- X_transcripts[X_transcripts %in% row.names(exp_matrix)]
A_transcripts_quantified <- A_transcripts[A_transcripts %in% row.names(exp_matrix)]

hist(log10(rowMeans(exp_matrix[A_transcripts_quantified,4:9])), col = 'green', ylab = 'log10 TPM', main = 'germ expression')
hist(log10(rowMeans(exp_matrix[X_transcripts_quantified,4:9])), add = T, col = 'red')
hist(log10(rowMeans(exp_matrix[L_transcripts_quantified,4:9])), add = T, col = 'yellow')

hist(log10(rowMeans(exp_matrix[A_transcripts_quantified,28:30])), col = 'green', ylab = 'log10 TPM', main = 'male embryos')
hist(log10(rowMeans(exp_matrix[X_transcripts_quantified,28:30])), add = T, col = 'red')
hist(log10(rowMeans(exp_matrix[L_transcripts_quantified,28:30])), add = T, col = 'yellow')

hist(log10(rowMeans(exp_matrix[A_transcripts_quantified,c(1:3,10:12)])), col = 'green', ylab = 'log10 TPM', main = 'no germ')
hist(log10(rowMeans(exp_matrix[X_transcripts_quantified,c(1:3,10:12)])), add = T, col = 'red')
hist(log10(rowMeans(exp_matrix[L_transcripts_quantified,c(1:3,10:12)])), add = T, col = 'yellow')



germ <- rowMeans(exp_matrix[L_transcripts_quantified,4:9])
no_germ <- rowMeans(exp_matrix[L_transcripts_quantified,c(1:3,10:12)])

length(germ)
length(no_germ)

plot(log10(germ) ~ log10(no_germ), pch = 20, cex = 0.2)

hist(log10(germ), col = 'yellow', probability = T)
hist(log10(no_germ), add = T, col = rgb(0.5,0.2,0.4, alpha = 0.3), probability = T)


mean(rowMeans(exp_matrix[L_transcripts_quantified,]) == 0)
sum(apply(exp_matrix[L_transcripts_quantified,4:9], 1, max) > 1)
sum(apply(exp_matrix[L_transcripts_quantified,4:9], 1, max) > 5)

expression_means <- apply(exp_matrix, 1, mean)
expression_max <- apply(exp_matrix, 1, max)

#png('figures/expression_max_logexpression.png')
#hist(log2(expression_max), breaks = 160, xlab = 'log2 expression', main = 'Sciara coprophila')
#dev.off()