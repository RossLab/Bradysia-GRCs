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

### Others

All the other species BUSCOs are here

```
../sciara_coprophila/21_L_age1/
```
