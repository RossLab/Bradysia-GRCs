## Gene content of L



### Hessian fly sex determination gene

The is a gene on a breakpoint of Hessian fly inversion that makes male-producing or female-producing females. The sequence is pasted bellow, and in a file `data/genome_wide_paralogy/Mdes_homologue_k05819.faa`. Rob found there is a copy of this gene on Sciara X, now we wonder if it's on L too.

```
QUERY=data/genome_wide_paralogy/Mdes_homologue_k05819.faa
GENOME=data/genome/genome.fasta
makeblastdb -in $GENOME -dbtype nucl

tblastn -query $QUERY -db $GENOME -evalue 1e-10 -outfmt 6 -num_threads 1 > data/genome_wide_paralogy/Mdes_homologue_k05819_2_Scop_genome.blast
```

Hits are: `NODE_3249`, `NODE_3029`, `NODE_3560`, `NODE_3527`, `NODE_8998`, using [this script](../scripts/L_gene_content/sex_determination_gene_k05819.R) we find that 3 of those are L, one is X and one is NA. However, the NA is quite well supported as L by kmers, only the coverage is quite strange, so it seems there are 3-4 L copies and one X copy. Also, the genes are fairly complete.

The sequence of the gene

```
>Mdes_homologue_k05819
MVVDSSKQNLQSENQKTNGNSKLGFSTEEALERLSSEIQDVLESYETKFKNDRLLSTVKSAFAKPYESPLSWFSILTTISVIVGLIVTKHYFTAICLLIILLMAVCLVIRESYLRHTEIHRKVRAAINEIELARQVCKDWTPENYPNVCSPMSPCVTLQLTYRDGKIVNLPWALLVNSDVIVMRPGSIAPTDCTELNGKIKFRAGEIYGLSQPIEAPVRPTIRHSLPDIICKLEKTPFLDNLVMTLDKAFKRPSTIYNQQRHLLITKYIQLWGVVICLLMTLLFAALHYNNKLYDKRLVELKWLNLFIELPLCSILPVIPMMFPVMWIIINLWGVARLQTYIEQPETIARPEQKRSFQEDLDTPTVECENVKLPFKEVFFNWIRLWNGDSHLLPRTANVVQVLGSVTALCCVDKKGILSWPNPTAEKVFFLRDSTEQTSNTGSIGSYESQERTSKSETCENSKKTCESSTQTPAYIDTKASITKRDSGRIVDTHKHPGTVAEVFDLTHDQCSPFRVDFDDHSWKKYIDSLKPLGLAILMNTCCQRTQEHYSKFCAHVTATAFFDKDLIPVTNRRCLCELAKQIGFTENAVKLFNLEGQIASYMHLQPDVVKRDIRFARQFQISSKVKVPFPHCLSVVMRDLQSNSLQMLTQGTADIVLDCCDDFWDGSDLRPLRPQERKRAQDFYQRNALTAYCTAFSYRPLRHGIAGALGETAYMELPPESTYKTSSHRDEYGGIFCTFDGTKSYIEPPNELQHSMSSDSLLFSDYKEDDVCDVDGCFEMQCHQVFIGMVTMQYQAQTDIVSLIERLERACIRFVHFSKENELRSRVFSEKMGLESGWNCHISLLSDEHSNAPSPQSSSRNFRAETCISNTKINIEDDIVVDSEMNQLLPPTSGLESSKNLSSSAPGAINEEFPSKSNKTEQSPNTSIESSQPKRFRVSKDSAMEEENCRSLSDYTDSTDQSAPVNFDMSNR
```