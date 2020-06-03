### making overview table

I link the Sciara BUSCO run

```
ln -s ../../sciara_coprophila/20_braker/busco.sciara.masked.spades.long/run_insecta_odb10/full_table.tsv data/busco.sciara.masked.spades.long_full_table.tsv
```

table of scaffold assignments

```
data/scaffold_assignment_tab_full.tsv
```

And now let's merge them running `Rscript scripts/make_BUSCO_scf_tab.R`. This will generate the table with assigments: [tables/BUSCO_assigned.tsv](https://github.com/RossLab/Sciara-L-chromosome/tables/BUSCO_assigned.tsv).

### Sorting the MAFFT alignments

For qtree the easiest way to specify an outgroup is to provide it as the first sequence in the fasta. The following script will require biopython (you can install that by pip) and resorts the fasta alignment so the outgroup is always the first (to standard output). If there is no outgroup no sequences are returned.

```
python3 scripts/Getting_L_age_using_BUSCO/sorting_MAFFT_alignmnets_of_BUSCO_genes.py <mafft.fasta> > <mafft_resorted.fasta>
```

### Others

All the other species BUSCOs are here

```
../sciara_coprophila/21_L_age1/
```
