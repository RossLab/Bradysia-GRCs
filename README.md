### Overview

*Bradysia coprophila* is a fungus gnat (fly) that carries germline restricted chromosomes (GRCs or 'L' chromosomes). We sequenced the germ tissue in *B. coprophila* and identified the GRCs using K-mer and coverage based techniques. We then compared the genes on the GRCs to those in the core genome (i.e. autosomes and X chromosome) to explore how the GRCs evolved in this species. The code in this repository was used in the analyses in this manuscript (https://doi.org/10.1101/2021.02.08.430288).

#### Note:
Within this repository, **L** is a synonym of **g**erm-line **r**estricted **c**hromosome (**GRC**).

### Project data

Later paths should be replaced by URLs in public reposutories.

- softmasked genome (Illumina): `/data/ross/mealybugs/analyses/sciara_coprophila/18_repeat/repeats/Scop_repeatmask/scop.spades2.min1kb.trimmed.fasta.masked`
- annotation: `/data/ross/mealybugs/analyses/sciara_coprophila/20_braker/braker.gff3`
- pacbio genome: `/data/ross/mealybugs/analyses/sciara_coprophila/Pacbio/4_racon2_sr/sciara.germ.pb.illumina1.rb.racon6pe.fasta`
- RNAseq of germ and soma (Males and Females):
`/data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/` and Mgerm/body or Fgerm/body directory
- Illumina raw reads:
`/data/ross/mealybugs/raw/transfer.genomics.ed.ac.uk/11372_Ross_Laura/raw_data/data_by_date/20180618/all_reads/`
- all vs all blasts
  - genes (nucleotides) `data/genome_wide_paralogy/genes_all_vs_all.blast`
  - proteins (transcribed genes) `data/genome_wide_paralogy/proteins_all_vs_all.blast`
- repeatmodeler/masker for genome annotation:
`/data/ross/mealybugs/analyses/Sciara-L-chromosome/data/repeats`

### Analysis documentation

Below is a list of manuscript sections containing analysis documentation and scripts used for those analyses

#### Genome assembly and annotation
  - associated data:
  * table linking Illumina contigs to annotated genes with assignments: `data/gene.assignment.tab.tsv`
    - associated r script: `scripts/gene.num.tab.R`
  * gene ID and mean coverage for that gene: `data/gene_cov_table.tsv`

#### Chromosome classification of the Sciara genome to GRC/X/A
  - [analyses documentation](analyses/assigment-of-L-X-A.md)
  - [associated scripts](scripts/kmer-assigment-of-L-X-A)
  - table with scores and assignments for all Illumina contigs
      `data/scaffold_assignment_tab_full.tsv`
      
#### Genome wide paralogy to identify GRC homologs
  - [analyses documentation](analyses/genome_wide_paralogy.md)


  - [scripts](scripts/paralog_divergence_prelim_analysis.R)
  - associated data
  * table with gene pairs in reciprocal blasts and gene cov, percent alignment between the blast pairs, length of genes, and assignments for genes. Before filtering based or gene/alignment length:  `data/ntgene_recip_blast_cov.tsv`
  * table with all info from paralog exploration, paralog ID, cov, type, and paralog freq after filtering based on gene/alignment length :  `data/filtered_paralog_tab.tsv`
    
#### Colinearity analysis
  - associated data
  * table of colinear blocks between chromosomes. Contains gene ID, chromosome assignment, scf in Illumina, scf in colinear block, block info, order in block and paralog ID, and assignment, and mean gene cov: `data/full_colinear_tab.tsv`
    - associated r script: `scripts/colinear_exploration.R`

#### GRC scaffold coverage analysis
  - [analyses documentation](analyses/cov_L_genes.md)


  - associated data
   * per gene coverage: `data/gene.cov.braker.annotation.tsv`

#### Phylogenetic placement of B. coprophila BUSCO genes
  - [analyses documentation](analyses/L_age_from_BUSCO.md)
  - [analyses documentation](analyses/Phylogeny_GRC_AX_paralogs_wBUSCO_IQtree.md)


  - associated scripts
   * [make_BUSCO_scf_tab.R](scripts/make_BUSCO_scf_tab.R)
   * [busco.sp.count.hist.R](scripts/busco.sp.count.hist.R)
  - associated data
   * [Table of assigment of BUSCO genes](tables/BUSCO_assigned.tsv)

#### GRC homolog identification in the _M. destructor_ (Cecidomyiidae) genome

  - [analyses documentation](analyses/GRC_homology_to_Mdestructor_vs_core_genome.md)


### Important places to look at

  - Google drive "Sciara_L_chromosome" (written documents; you should have access if you are supposed to have an access)
  - [Code](https://github.com/RossLab/Sciara-L-chromosome/blob/master/scripts/kmer-assigment-of-L-X-A/L-assignment.R#L28-L30) that **decides about assigments** what is L / Lc / Lk / Lp / AX / AXp (we can latter one separate A and X using either male coverage or the reference genome). In total the assigment fractions now are (Mbp):

  |    L    |    Lc    |    A    |    Ac    |    X    |    Xc    |    NA    |
  | ------- | -------- | ------- | -------- | ------- | -------- | -------- |
  | 154.1   |  6.8     |  162.4  |  9       |  52.9   |  2.4     |  10.1    |

