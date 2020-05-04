
Here should be some description in future.

#### Notes

Within this repository, **L** is a synonym of **g**erm-line **r**estricted **c**hromosome (**GRC**).

### Data

Later paths should be replaced by URLs in public reposutories.

- softmasked genome (Illumina): `/data/ross/mealybugs/analyses/sciara_coprophila/18_repeat/repeats/Scop_repeatmask/scop.spades2.min1kb.trimmed.fasta.masked`
- annotation: `/data/ross/mealybugs/analyses/sciara_coprophila/20_braker/braker.gff3`
- pacbio genome: `/data/ross/mealybugs/analyses/sciara_coprophila/Pacbio/4_racon2_sr/sciara.germ.pb.illumina1.rb.racon6pe.fasta`
- RNAseq of germ and soma (Males and Females):
`/data/ross/mealybugs/raw/11791_Ross_Laura/raw_data/20190812/` and Mgerm/body or Fgerm/body directory
- Illumina raw reads:
`/data/ross/mealybugs/raw/transfer.genomics.ed.ac.uk/11372_Ross_Laura/raw_data/data_by_date/20180618/all_reads/`


### Subprojects

This in an incomplete list, only those there are here on GitHub

#### Sorting out the Sciara genome to L/X/A using coverage differences between germ and somatic libs
  - [analyses documentation](scripts/coverage_assignment_notes.md)

#### Sorting out the Sciara genome to L/X/A using kmers
  - [project board](https://github.com/orgs/RossLab/projects/1)
  - [analyses documentation](analyses/kmer-assigment-of-L-X-A.md)
  - [associated scripts](scripts/kmer-assigment-of-L-X-A)

#### Using Busco genes to determine the age of the L chromosome

#### Genome wide paralogy
  - [analyses documentation](analyses/genome_wide_paralogy.md)
  
#### Expression of L chromosome genes


### Important places to look at

  - Google drive "Sciara_L_chromosome" (written documents; you should have access if you are supposed to have an access)
  - [Code](https://github.com/RossLab/Sciara-L-chromosome/blob/master/scripts/kmer-assigment-of-L-X-A/L-assignment.R#L28-L30) that **decides about assigments** what is L / Lc / Lk / Lp / AX / AXp (we can latter one separate A and X using either male coverage or the reference genome). In total the assigment fractions now are (Mbp):

  |    L    |    Lc    |    A    |    Ac    |    X    |    Xc    |    NA    |
  | ------- | -------- | ------- | -------- | ------- | -------- | -------- |
  | 154.1   |  6.8     |  162.4  |  9       |  52.9   |  2.4     |  10.1    |

  - [Table of assigment of BUSCO genes](tables/BUSCO_assigned.tsv)
