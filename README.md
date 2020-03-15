
Here should be some description in future.

Notes:
  - Within the repository, **L** is a synonim of **GRC**

### Projects

#### Sorting out the Sciara genome to L/X/A using kmers:
  - [project board](https://github.com/orgs/RossLab/projects/1)
  - [analyses documentation](analyses/kmer-assigment-of-L-X-A.md)
  - [associated scripts](scripts/kmer-assigment-of-L-X-A)

### Important places to look at

  - [Code](https://github.com/RossLab/Sciara-L-chromosome/blob/master/scripts/kmer-assigment-of-L-X-A/L-assignment.R#L28-L30) that **decides about assigments** what is L / Lc / Lk / Lp / AX / AXp (we can latter one separate A and X using either male coverage or the reference genome). In total the assigment fractions now are (Mbp):

  |    L    |    Lk    |    Lc    |    Lp    |    AX    |    AXp    |
  | ------- | -------- | -------- | -------- | -------- | --------- |
  | 129.4   |  30.9    |  4.4     |  6       |  217.6   |  9.4      |


  - [Table of assigment of BUSCO genes](tables/BUSCO_assigned.tsv)
