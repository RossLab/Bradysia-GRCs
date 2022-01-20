


### Objective
We want to explore the origin of GRC genes by looking at their phylogenetic placement in phylogenies made with BUSCO genes from representitive species spanning Sciaroidea (the group that contains _B. coprophila_ and other Sciaridae, Cecidomyiidae (another fly group that has GRCs) and other gnats. We're going to use BUSCO genes for this analysis since we can easily identify homologs in all the species we're interested in.

### Analysis
#### Before step 1:
- Need to add some additional detail about how we ran BUSCO on the gnat genomes and pulled out homologs based on BUSCO ID.
- the BUSCO IDs in the 4b_busco_allsp_aa2/ folder are only the ID's with more than 80 % of the species in the alignment. .

#### Counting how many species are in BUSCO ID:
```
for gene in $(cat s.coprophila.duplicated.ID.txt); do
  count=$(grep ">" 4_busco_allsp_aa2/$gene | cut -f 1,2 -d "_" | sort | uniq | wc -l)
  echo -e "$gene \t $count" >> sp.count.dup.busco.txt
  done
for gene in $(cat s.coprophila.single.ID.txt); do
  count=$(grep ">" 4b_busco_allsp_aa2/$gene | cut -f 1,2 -d "_" | sort | uniq | wc -l)
  total=$(grep ">" 4b_busco_allsp_aa2/$gene | cut -f 1,2 -d "_" | wc -l)
  echo -e "$gene \t $count \t $total">> sp.count.sc.busco.uniq.tot.txt
  done
```

```
qsub -o logs -e logs -cwd -N spcount -V -pe smp64 1 -b yes './command.sp.count.sh'
```
Now I'm going to import these tables into R and make some histograms of which aa seqs to not include

rscript:
`/data/ross/mealybugs/analyses/Sciara-L-chromosome/scripts/busco.sp.count.hist.R`
If we include all of the genes that have more than 80% species in them then we include anything over 12.8 sp which turns out to be 1200 genes (I think that's good).

busco.gene.80.txt (genes we want to keep)
```
for gene in $(cat busco.gene.80.txt); do
  mv 4_busco_allsp_aa2/$gene 4b_busco_allsp_aa2/$gene
done
```
#### Moving genes of interest to one folder
- IDs that have 80% of species present
```
for faa_file in $(cat data/phylogenies_grc/BUSCO_80_all_id.txt); do
  mv /data/ross/mealybugs/analyses/sciara_coprophila/21_L_age1/4b_busco_allsp_aa2/"$faa_file".faa data/phylogenies_grc/other_aa_seq_allsp/"$faa_file".faa
done
```
#### Keeping only longest orf in BUSCO amino acid fastas

location: `/data/ross/mealybugs/analyses/sciara_coprophila/21_L_age1/4b_busco_allsp_aa2/`
script: `scripts/get_lengest_busco_sp.py`

Script to run all BUSCO aa fasta files through longest orf script, which pulls out the longest orf for species other than B. cop and labels B. cop so that the A/X copy is labelled A and the GRC copy or copies are labelled L
```
for faa_file in data/phylogenies_grc/other_aa_seq_allsp/aa_seq/*; do
  gene=$(echo "$faa_file" | cut -f 5 -d "/" | cut -f 1 -d ".")
  python3 /data/ross/mealybugs/analyses/Sciara-L-chromosome/scripts/get_lengest_busco_sp.py "$faa_file" > data/phylogenies_grc/other_aa_seq_allsp/aa_seqs_longestorf/"$gene".fasta
done
```

#### Running mafft on each fasta (BUSCO ID) to get alignment

##### Running mafft on each fasta to get alignment
```
for faa_file in data/phylogenies_grc/other_aa_seq_allsp/aa_seqs_longestorf/*.fasta; do
  gene=$(echo "$faa_file" | cut -f 5 -d "/" | cut -f 1 -d ".")
  linsi "$faa_file" > data/phylogenies_grc/other_aa_seq_allsp/mafft_alignments/"$gene".fasta
done
```

#### Reorganising mafft alignment so outgroup is first species
*outgroup has to be first species in alignment*
script to reformat mafft files
`/data/ross/mealybugs/analyses/Sciara-L-chromosome/scripts/Getting_L_age_using_BUSCO/sorting_MAFFT_alignmnets_of_BUSCO_genes.py`

```
for fasta_file in data/phylogenies_grc/other_aa_seq_allsp/mafft_alignments/*.fasta; do
  gene=$(echo "$fasta_file" | cut -f 5 -d "/" | cut -f 1 -d ".")
  python3 /data/ross/mealybugs/analyses/Sciara-L-chromosome/scripts/Getting_L_age_using_BUSCO/sorting_MAFFT_alignmnets_of_BUSCO_genes.py "$fasta_file" > data/phylogenies_grc/other_aa_seq_allsp/mafft_sorted/"$gene"_sorted.fasta
done
```
#### Making gene trees with IQtree for each BUSCO id
this is making a phylogeny for each BUSCO ID separately
```
for fasta_file in data/phylogenies_grc/other_aa_seq_allsp/mafft_sorted/*; do
  gene=$(echo "$fasta_file" | cut -f 5 -d "/" | cut -f 1 -d "_")
  iqtree -s "$fasta_file" -alrt 1000 -bb 1000 -nt AUTO -ntmax 8 -pre data/phylogenies_grc/other_aa_seq_allsp/genetrees_iqtree/"$gene"
done
```
- We manually inspected some of the genetrees and it seems that most GRC genes either are in the Sciaridae clade in the phylogeny or in the Cecidomyiidae clade. We made a script to summarise which of these two clades (or other) each _B. coprophila_ BUSCO gene is placed in. It also gives us some extra information like the bootstrap value at the nearest node for _B. coprophila_ genes and the branch length of the _B. coprophila_ branches.

Script: `scripts/phylogeny/treefile2table_of_neighbors.py`

```
qsub -o logs -e logs -cwd -N phy_summary -V -pe smp64 1 -b yes 'python3 scripts/phylogeny/treefile2table_of_neighbors.py data/phylogenies_grc/other_aa_seq_allsp/genetrees_iqtree/ > tables/phylogeny_summary_allotherBUSCOS.tsv'
```

Update 7.1.2022. The table now contains `scfs` column as well, which will specify the original scaffold where each of the BUSCO copies was placed. Like this the table can be regenerated

```
python3 scripts/phylogeny/treefile2table_of_neighbors.py data/phylogenies_grc/all_buscos/genetrees_iqtree > tables/phylogeny_summary_allBUSCOS_redone.tsv
```

- next, summarising the output of this script in R:
script: `scripts/BUSCO_all_plots.R`

#### Concatenated phylogenies for GRC-core gene BUSCO genes:
I'm going to make two different concatentated phylogenies. The first will be for the genes that fall within the Cecidomyiidae and the second will be for the ones that fall within the Sciaridae. Then I'll put an N= on the figure showing how many genes went into each. I'm going to change the order of figure 5 so it's the gene trees first, then the concatentated trees.

##### Extracting GRC-Core gene BUSCOs and separating based on phylogenetic placement. 
- Then taking away scaffold info from species name so when I run IQtree genes from the same species will concatenate
```
for faa_file in $(cat data/phylogenies_grc/all_buscos/ceci_buscos_AL_ids.tsv); do
  cp data/phylogenies_grc/all_buscos/mafft_sorted/"$faa_file"_sorted.fasta data/phylogenies_grc/all_buscos/ceci_busco_tree/"$faa_file"_sorted.fasta
done
for faa_file in $(cat data/phylogenies_grc/all_buscos/sciaridae_buscos_AL_ids.tsv); do
  cp data/phylogenies_grc/all_buscos/mafft_sorted/"$faa_file"_sorted.fasta data/phylogenies_grc/all_buscos/sciaridae_busco_tree/"$faa_file"_sorted.fasta
done
```
```
for faa_file in *; do
	cat "$faa_file" | awk 'BEGIN { FS = "_" }; />/{ print $1"_"$2 } !/^>/ { print $0 }' > renamed."$faa_file" && mv renamed."$faa_file" "$faa_file"
done
```
#### Making the phylogenies
```
qsub -o logs -e logs -cwd -N concat_phy -l h=biggar -V -pe smp64 32 -b yes 'iqtree -s data/phylogenies_grc/all_buscos/sciaridae_busco_tree/ -alrt 1000 -bb 1000 -nt AUTO -ntmax 32 -mset LG -msub nuclear -pre data/phylogenies_grc/all_buscos/concat_phys/sciaridae_AL_tree'
qsub -o logs -e logs -cwd -N concat_phy -l h=bigbang -V -pe smp64 32 -b yes 'iqtree -s data/phylogenies_grc/all_buscos/ceci_busco_tree/ -alrt 1000 -bb 1000 -nt AUTO -ntmax 32 -mset LG -msub nuclear -pre data/phylogenies_grc/all_buscos/concat_phys/ceci_AL_tree'
```


#### Amino acid sequence composition of genes

Looking at whether weird amino acid composition of GRC genes could account for phylogenetic position (due to long branch attraction).

`data/phylogeny_mafft`

```python
import itertools
from collections import defaultdict
from os import listdir
from os.path import isfile, join
from sys import stdout

sp2aas = defaultdict(lambda: defaultdict(int))

align_dir = 'data/phylogenies_grc/all_buscos/mafft_sorted/'
filenames = [join(align_dir, f) for f in listdir(align_dir) if isfile(join(align_dir, f))]

for filename in filenames:
  with open(filename) as file:
    for header,sequence in itertools.zip_longest(*[file]*2):
      species = "_".join(header[1:].split('_')[0:2])
      sequence = sequence.rstrip('\n').replace("-", "")
      for aa in sequence:
        sp2aas[species][aa] += 1

# nt_aln_dir = 'data/phylogeny_mafft'

all_aas = set.intersection(*map(set,[sp2aas[sp].keys() for sp in sp2aas.keys()]))
stdout.write("Species\t" + '\t'.join(all_aas) + '\n')

for sp in sp2aas.keys():
  aa_counts = [sp2aas[sp][aa] for aa in all_aas]
  aa_freq = [str(round(aa_c / sum(aa_counts), 4)) for aa_c in aa_counts]
  stdout.write(sp + "\t" + '\t'.join(aa_freq) + '\n')
```

I generated `tables/BUSCO_aa_composition.tsv` table (`/data/ross/mealybugs/analyses/Sciara-L-chromosome/tables/BUSCO_aa_composition.tsv`) with aa relative composition and now I can plot similarities between samples.

```R
aa_comp <- read.table('tables/BUSCO_aa_composition.tsv', header = T)

# The mtcars dataset:
data <- as.matrix(aa_comp[, -1])
row.names(data) <- aa_comp[, 1]

# Default Heatmap
heatmap(data)
```

