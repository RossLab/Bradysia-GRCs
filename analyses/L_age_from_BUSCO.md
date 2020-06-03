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
