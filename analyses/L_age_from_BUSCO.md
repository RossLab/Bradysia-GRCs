What we have at the moment (I'll add the appropriate code later on at some point)

1. there is a large number of BUSCO genes that are either on the L chromosome or have a paralog on the L chromosome (i.e. in some way associated with L).This is interesting for two reasons:
- it seems like there are a lot of duplicated BUSCO genes on the L chromosome
- there are some single copy L BUSCO's. this is unexpected if there are no autosomal paralogs since the BUSCO genes are conserved genes and we wouldn't expect them to be on chromosomes only present in germ tissue/ with a weird chromosome cycle
2. if you put the BUSCO genes in a phylogeny you see that the L gene copies are either closely related to the Hessian fly (Mayetolia destructor), or is closely related to S. coprophila.
3. We have some alignments and phylogenies for all BUSCO genes for S. coprophila and also for all other Sciaroidea species we have genomes for (from Noelle's paper). We need to curate these alignments (make some rules about how many species have to be present and maybe trim the alignments at places where fewer than a certain number of species have sequence. Then we want to re-run phylogenies with an outgroup in the phylogeny

## Status:
- I've decided to take the species exceria fusca out of the phylogenies because it has a 50% or so BUSCO score. I'm guessing the assembly isn't that good

directory I'm runing code from:
`/data/ross/mealybugs/analyses/sciara_coprophila/21_L_age1/`
directory to put aa seq into:
`4_busco_allsp_aa2/`
list of all busco aa id's from s.cop: s.coprophila.allaa.ID.txt

### Moving all busco aa seqs (that have complete copies in S.cop) and all other sciaridae sp (except e. fusca which has only a 50% busco score so we've excluded it) to one directory according to gene:
```
for gene in $(cat s.coprophila.allaa.ID.txt); do
	for faa_file in */run_insecta_odb10/busco_sequences/*_copy_busco_sequences/"$gene"; do
		species=$(echo "$faa_file" | cut -f 1 -d "/" | cut -f 2 -d "." )
		cat "$faa_file" | awk -v sp="$species" 'BEGIN { FS = "\t" }; />/{ print ">" sp "_" substr($1,2) } !/^>/ { print $0 }' >> 4_busco_allsp_aa2/"$gene"
	done
done
```
```
qsub -o logs -e logs -cwd -N sort -V -pe smp64 2 -b yes './command.busco.allsp.sh'
```
### Counting how many species are in each gene list:
```
for gene in $(cat s.coprophila.duplicated.ID.txt); do
  count=$(grep ">" 4_busco_allsp_aa2/$gene | cut -f 1,2 -d "_" | sort | uniq | wc -l)
  echo -e "$gene \t $count" >> sp.count.dup.busco.txt
  done
for gene in $(cat s.coprophila.single.ID.txt); do
  count=$(grep ">" 4_busco_allsp_aa2/$gene | cut -f 1,2 -d "_" | sort | uniq | wc -l)
  echo -e "$gene \t $count">> sp.count.sc.busco.txt
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
### Aligning all gene seqs in 4b_busco_allsp_aa2/ (using mafft linsi)
```
for faa_file in 4b_busco_allsp_aa2/*.faa; do
  gene=$(echo "$faa_file" | cut -f 2 -d "/" | cut -f 1 -d ".")
  linsi "$faa_file" > 5_mafft2/"$gene".fasta
done
```
```
qsub -o logs -e logs -cwd -N mafft -V -pe smp64 1 -b yes './command.mafft_3.sh'
```
ok, done. I've decided not to do any additional filtering of alignments.

### now need to go through and figure out what outgroup we should use for each tree. if possible, use sylvicola_fuscatus, but penthetria_funebris is also an outgroup, can I use both for the outgroup and what is the easiest way to tally which ones have which outgroups and then set the iqtree command

-o is the way to specify the outgroup, I think we'll need to rename the nodes to just the species names to use this those since I think the names have to be exact (I can try and see though)

just testing to see how exact - o name has to be
qsub -o logs -e logs -cwd -N iqtree -V -pe smp64 8 -b yes 'iqtree -s 5_mafft2/147143at50557.fasta -alrt 1000 -bb 1000 -nt AUTO -ntmax 8 -pre 147143at50557 -o sylvicola_fuscatus'
ok, name needs to be exact so I will need to go through and take off all bits after species name in each fasta file.


below is some code that doesn't work. We have a problem that iqtree isn't recognizing the -o outgroup name. It keeps saying that it the species isn't in the list even though it definitely is...


e.g name: >gnoriste_bilineata_k141_243800:598-1264 <unknown description>

```
for fasta_file in 5_mafft2/*; do
  gene=$(echo "$fasta_file" | cut -f 2 -d "/" | cut -f 1 -d ".")
  outgroup=$(cat "$fasta_file" | grep '>sylvicola_fuscatus' | cut -f2 -d ">" | cut -f1 -d " ")
  [[ -z "$outgroup" ]] && outgroup=$(cat "$fasta_file" | grep '>penthetria_funebris' | cut -f2,3 -d ">")
  iqtree -s "$fasta_file" -alrt 1000 -bb 1000 -nt AUTO -ntmax 8 -pre "$gene" -o "$outgroup"
done

for fasta_file in 5_mafft2/73852at50557.fasta; do
  gene=$(echo "$fasta_file" | cut -f 2 -d "/" | cut -f 1 -d ".")
  outgroup=$(cat "$fasta_file" | grep '>sylvicola_fuscatus')
  [[ -z "$outgroup" ]] && outgroup=$(cat "$fasta_file" | grep '>penthetria_funebris')
  # deletes the > at the beginning of the fasta entry name
  outgroup=${outgroup#?};
  iqtree -s "$fasta_file" -alrt 1000 -bb 1000 -nt AUTO -ntmax 8 -pre "$gene" -o "$outgroup"
done

qsub -o logs -e logs -cwd -N iqtree -V -pe smp64 8 -b yes './command.iqt

qsub -o logs -e logs -cwd -N iqtree -V -pe smp64 8 -b yes 'iqtree -s 5_mafft2/73852at50557_2.fasta -o trichosia_splendens_k141_213139 -alrt 1000 -bb 1000 -nt AUTO -ntmax 8 -pre 73852at50557_2.fasta'
```
try above before moving on...

```
awk '/^>/ {$0=$1} 1' file.fasta

sed '/^>/ s/ .*//' 5_mafft2/73852at50557.fasta > 5_mafft2/73852at50557_2.fasta

for gene in 5_mafft2/*; do
sed '/^>/ s/ .*//' "$gene" > 5_mafft2/"$gene"
done


100236at50557.faa

>sylvicola_fuscatus_k141_82067:1-4813
>penthetria_funebris_k141_62309:15-4604
>bolitophila_cinerea_k141_24180:6608-9217
>bolitophila_hybrida_k141_22529:3316-11503
>catotricha_subobsoleta_k141_662152:1-8407
>diadocidia_ferruginosa_k141_263832:1-7489
>gnoriste_bilineata_k141_310207:626-4918
>lestremia_cinerea_k141_2488:1-5748
>macrocera_vittata_k141_140786:5348-8357
>mayetiola_destructor_GL503022.1:66351-71763
>phytosciara_flavipes_k141_241131:2482-8665
>platyura_marginata_k141_36182:4983-13970
>sciara_coprophila_NODE_747_length_56249_cov_136.969967:16800-31500
>symmerus_nobilis_k141_6397:1726-7451
>trichosia_splendens_k141_426474:1-5856
>porricondyla_nigripennis_
```

### Previous code

I am not sure now where this code should be placed. The first bit seems to belong to point 1. in the overview (exploratory analysis of BUSCO genes). The second part is for sorting MAFFT alignments for phylogeny, not sure if it was used at all, so I keep it here to be sorted out later:

#### Overview of BUSCO genes table

I link the Sciara BUSCO run

```
ln -s ../../sciara_coprophila/20_braker/busco.sciara.masked.spades.long/run_insecta_odb10/full_table.tsv data/busco.sciara.masked.spades.long_full_table.tsv
```

table of scaffold assignments

```
data/scaffold_assignment_tab_full.tsv
```

And now let's merge them running `Rscript scripts/make_BUSCO_scf_tab.R`. This will generate the table with assigments: [tables/BUSCO_assigned.tsv](https://github.com/RossLab/Sciara-L-chromosome/tables/BUSCO_assigned.tsv).

#### Sorting the MAFFT alignments

For qtree the easiest way to specify an outgroup is to provide it as the first sequence in the fasta. The following script will require biopython (you can install that by pip) and resorts the fasta alignment so the outgroup is always the first (to standard output). If there is no outgroup no sequences are returned.

```
python3 scripts/Getting_L_age_using_BUSCO/sorting_MAFFT_alignmnets_of_BUSCO_genes.py <mafft.fasta> > <mafft_resorted.fasta>
```

### FigTree Notes


```bash
phylopath=/Volumes/dump/projects/PGE/Sciara-L-Chromosome/data/GRC_phylogenies
for treefile in $phylopath/15_genetrees_L/*treefile; do
  pngfile=$(basename $treefile .treefile).png
  java -Xms64m -Xmx512m -jar lib/figtree.jar -graphic PNG $treefile $phylopath/plot_$pngfile
done
```







*Christina's note: Code below is to run Iqtree on all BUSCO id's with one GRC gene and one A/X gene, the code below was what worked in the end, I think some of the code above wasn't used*


note: had some issues initially so started again in 10_seqs2. need to delete 1_seq_etc and 2_ 3_ any others before 10_

note2: realised that we can't do species tree with 2L genes since we don't know which L each gene belongs to so we can only look at gene trees for those ones (not that many anyways)


**Objective of script** is to make a species alignment and phylogeny for A/X and L duplicate buscos (there are more than 300) and gene trees for each gene, then to try to reconcile the species tree and gene trees with concordance factors in IQtree.


### Before step 1: 
- the BUSCO IDs in the 4b_busco_allsp_aa2/ folder are only the ID's with more than 80 % of the species in the alignment. Might need to take out any without all the species later.

### Step1: Keeping only longest orf in BUSCO aa fastas

location: `/data/ross/mealybugs/analyses/sciara_coprophila/21_L_age1/4b_busco_allsp_aa2/`
script: `scripts/get_lengest_busco_sp.py`

Script to run all BUSCO aa fasta files through longest orf script, which pulls out the longest orf for species other than S. cop and labels S. cop so that the A/X copy is labelled A and the L copy/ies are labelled L (need to run in channotation env, chodson python3 doesnt work)
```
for faa_file in 4b_busco_allsp_aa2/*.faa; do
  gene=$(echo "$faa_file" | cut -f 2 -d "/" | cut -f 1 -d ".")
  python3 /data/ross/mealybugs/analyses/Sciara-L-chromosome/scripts/get_lengest_busco_sp.py "$faa_file" > 7_GRC_phylogenies/1_seqs_longestorfs/"$gene".fasta
done
```
for faa_file in 4b_busco_allsp_aa2/10008at50557.faa; do
  gene=$(echo "$faa_file" | cut -f 2 -d "/" | cut -f 1 -d ".")
  python3 /data/ross/mealybugs/analyses/Sciara-L-chromosome/scripts/get_lengest_busco_sp.py "$faa_file" > 7_GRC_phylogenies/"$gene".fasta
done

```
qsub -o logs -e logs -cwd -N longorf -V -pe smp64 1 -b yes './command.longestorf.sh'
```

### Step2: running mafft on each fasta to get alignment

```
for faa_file in 7_GRC_phylogenies/10_seqs2/*.fasta; do
  gene=$(echo "$faa_file" | cut -f 3 -d "/" | cut -f 1 -d ".")
  linsi "$faa_file" > 7_GRC_phylogenies/2_mafft_alignments/"$gene".fasta
done
```
```
qsub -o logs -e logs -cwd -N mafft -V -pe smp64 1 -b yes './command.mafft_grc.sh'
```

### Step3: reorganising mafft alignment so outgroup is first species
*outgroup has to be first species in alignment*
script to reformat mafft files
`/data/ross/mealybugs/analyses/Sciara-L-chromosome/scripts/Getting_L_age_using_BUSCO/sorting_MAFFT_alignmnets_of_BUSCO_genes.py`

```
for fasta_file in 7_GRC_phylogenies/12_mafft_alignments/*.fasta; do
  gene=$(echo "$fasta_file" | cut -f 3 -d "/" | cut -f 1 -d ".")
  python3 /data/ross/mealybugs/analyses/Sciara-L-chromosome/scripts/Getting_L_age_using_BUSCO/sorting_MAFFT_alignmnets_of_BUSCO_genes.py "$fasta_file" > 7_GRC_phylogenies/13_mafft_sorted/"$gene"_sorted.fasta
done
```
```
qsub -o logs -e logs -cwd -N sorting -V -pe smp64 2 -b yes './command.mafft.sort.sh'
```
here

### Step4: take out BUSCOs with one GRC gene and one A/X
took out BUSCO id's for BUSCO groups:
- A-L and L-X: `data/Lother_paralogs.txt`
- A-L-L and L-L-X: `data/LLother_paralogs.txt`

script to do this in R (probably should be separate script): `make_BUSCO_scf_tab`
then 
```
cut -f1 Lother_paralogs.txt | uniq > L_XA_buscoid.Illumina.txt
```

*note* need to redo file that gives BUSCO gene pairs with Illumina BUSCO file
```
for gene in $(cat L_XA_buscoid.Illumina.txt); do
  cp 7_GRC_phylogenies/13_mafft_sorted/"$gene"_sorted.fasta 7_GRC_phylogenies/14_mafft_L_other/
done
for gene in $(cat LL_XA_buscoid.Illumina.txt); do
  cp 7_GRC_phylogenies/13_mafft_sorted/"$gene"_sorted.fasta 7_GRC_phylogenies/14b_mafft_LL_other/
done
```
output (BUSCO's not present in 80% of species in phylogeny, we omitted those BUSCO ids)
```
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/2418at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/12466at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/13582at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/18400at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/46968at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/47520at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/57218at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/59099at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/59985at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/61099at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/66657at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/68581at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/70521at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/76908at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/77723at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/78674at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/81919at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/84147at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/87518at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/89135at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/91691at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/97917at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/102539at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/106023at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/106197at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/106896at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/111975at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/119324at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/123925at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/128559at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/137014at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/146376at50557_sorted.fasta': No such file or directory
```
(32)
```
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/39740at50557_sorted.fasta': No such file or directory
cp: cannot stat '7_GRC_phylogenies/13_mafft_sorted/93865at50557_sorted.fasta': No such file or directory
```
(2)
```
qsub -o logs -e logs -cwd -N sorting -V -pe smp64 2 -b yes './command.mvalignments.sh'
```

### Step5: running IQtree on all alignments to make species tree and also on all BuscoID's separately to make gene trees.
This is renamming the species names in the fasta file to make the species tree. I think I'll do this on all the files in the 14_ directory since they're copied from 13_ anyways
#### For species trees
```
for faa_file in *; do
	cat "$faa_file" | awk 'BEGIN { FS = "_" }; />/{ print $1"_"$2 } !/^>/ { print $0 }' > renamed."$faa_file" && mv renamed."$faa_file" "$faa_file"
done
```
```
qsub -o logs -e logs -cwd -N iqtree -V -pe smp64 8 -b yes 'iqtree -s 7_GRC_phylogenies/14_mafft_L_other/ -alrt 1000 -bb 1000 -nt AUTO -ntmax 8 -pre 7_GRC_phylogenies/15c_sp_tree_Lother/Lother_sp'
```

For species trees (in this case species tree includes the GRC(s) as a "species")


#### For gene trees
```
for fasta_file in 7_GRC_phylogenies/14_mafft_L_other/*; do
  gene=$(echo "$fasta_file" | cut -f 3 -d "/" | cut -f 1 -d "_")
  iqtree -s "$fasta_file" -alrt 1000 -bb 1000 -nt AUTO -ntmax 8 -pre 7_GRC_phylogenies/15_genetrees_L/"$gene"
done
for fasta_file in 7_GRC_phylogenies/14b_mafft_LL_other/*; do
  gene=$(echo "$fasta_file" | cut -f 3 -d "/" | cut -f 1 -d "_")
  iqtree -s "$fasta_file" -alrt 1000 -bb 1000 -nt AUTO -ntmax 8 -pre 7_GRC_phylogenies/15b_genetrees_LL/"$gene"
done
```

```
qsub -o logs -e logs -cwd -N iqtree -V -pe smp64 8 -b yes './command.iqtree_L.sh'
```

### to compute concordance factors
command from iqtree site
```
iqtree -t concat.treefile --gcf loci.treefile -p ALN_DIR --scf 100 --prefix concord -T 10
```
```
qsub -o logs -e logs -cwd -N iqtree -V -pe smp64 10 -b yes 'iqtree -t 7_GRC_phylogenies/15c_sp_tree_Lother/Lother_sp.treefile --gcf 7_GRC_phylogenies/15_genetrees_L/*.treefile -p 7_GRC_phylogenies/14_mafft_L_other/ --scf 100 --prefix 7_GRC_phylogenies/16_concord/concord -T 10'
```
