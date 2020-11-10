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
output
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
```
for faa_file in *; do
	cat "$faa_file" | awk 'BEGIN { FS = "_" }; />/{ print $1"_"$2 } !/^>/ { print $0 }' > renamed."$faa_file" && mv renamed."$faa_file" "$faa_file"
done
```
```
qsub -o logs -e logs -cwd -N iqtree -V -pe smp64 8 -b yes 'iqtree -s 7_GRC_phylogenies/14_mafft_L_other/ -alrt 1000 -bb 1000 -nt AUTO -ntmax 8 -pre 7_GRC_phylogenies/15c_sp_tree_Lother/Lother_sp'
```

For species trees (in this case species tree includes the GRC(s) as a "species")


For gene trees
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

### Now need to compute concordance factors
command from iqtree site
```
iqtree -t concat.treefile --gcf loci.treefile -p ALN_DIR --scf 100 --prefix concord -T 10
```
```
qsub -o logs -e logs -cwd -N iqtree -V -pe smp64 10 -b yes 'iqtree -t 7_GRC_phylogenies/15c_sp_tree_Lother/Lother_sp.treefile --gcf 7_GRC_phylogenies/15_genetrees_L/*.treefile -p 7_GRC_phylogenies/14_mafft_L_other/ --scf 100 --prefix 7_GRC_phylogenies/16_concord/concord -T 10'
```