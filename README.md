# Hornwort Genomes Paper

## Methods
1. Sampling and Cultivation
2. DNA and RNA Sequencing
3. Genome Assembly and Scaffolding
4. Annotation
	a. Repetitive Elements
	b. Gene Models
	c, Methylation
5. Chromosome Structure
6. Orthogroup Inference
7. Gene Family Expansion/Contraction
8. Synteny
	a. Within hornworts
	b. Within bryophytes
	c. Chromosome evolution
9. Pangenome
10. WGD
11. Horizontal Gene Transfer
12. Stomatal Gene Expression
13. Cyanobacterial Symbiosis
14. GID1-DELLA
15. Anthocyanin pathway
16. Pyrenoid genes
17. LTR taxonomy, distributions, ages

### Sampling and Cultivation
Sporophytes were collected from the wild (Table 1). Depending on size, one to three sporophytes were macerated with forceps in 1.5 mL centrifuge tubes containing 1 mL of 0.1% Tween in sterile water, vortexed for 30 seconds, then centrifuged at 5000 RPM for 5 minutes to pull down spores. Sporophyte debris was removed, then the supernatant pipetted off. One mL of sterile water was added to the spore pellet, vortexed for 30 seconds, then centrifuged as above and the water removed. This washing was repeated three times. Spores were resuspended in 30 uL of sterile water and split into 10 uL aliquots. Aliquots were added to 50 uL of 2.5%, 5%, and 10% dilutions of commercial bleach, vortexed for five seconds, and let sit for one minute. After one minute, 50 uL of sterile 0.1 M sodium thiosulfate was added and immediately vortexed for five seconds to quench the reaction.

| Taxon | Family | ID | Locality |
|-------|--------|----|----------|
| Anthoceros agrestis 'Bonn' | Anthocerotaceae | Bonn | Bonn, Switzerland |
| A. agrestis 'Oxford' | Anthocerotaceae | Oxford | Oxford, UK? |
| A. fusiformis | Anthocerotaceae | Chico A |  |
| A. punctatus | Anthocerotaceae | Meeks |  |
| Leiosporoceros dussii | Leiosporocerotaceae | K | Panama? |
| Megaceros flagellaris | Dendrocerotaceae | A |  Hawaii |
| Notothylas orbicularis | Notothyladaceae |  |   |
| Paraphymatoceros pearsonii | Notothyladaceae |   |
| Phaeoceros carolinianus | Notothyladaceae |  |   |
| Phaeoceros sp. | Notothyladaceae |  |   |
| Phymatoceros phymatodes | Phymatocerotaceae |  |  California |

### DNA and RNA Sequencing

### Genome Assembly and Scaffolding
ONT reads less than 5 kbp were removed and the remainder were assembled with Flye v2.9. Contigs were corrected with Illumina DNA sequence data using Pilon v1.24. HiC libraries were prepared, sequenced, and scaffolded by Phase Genomics (Seattle, WA). TGS-Gapcloser was used to fill gaps between scaffolds with ONT reads and polish filled gaps with Illumina data.

### Annotation
#### Repetitive Elements
Custom repeat libraries were constructed for each species using EDTA v2. The TE library output by EDTA was filtered by extracting LTRs labeled as unknown and BLASTing their nucleotide sequences against a database of transposases. Any sequences with significant hits were retained as part of the TE library, and the others removed. The TE library was then BLASTed against known plant proteins in Uniprot  


#### Gene Models


#### Methylation


### Orthogroups
Orthogroups inferred at three taxonomic levels: Anthocerotaceae, Byrophyta, and Viridiplantae. Hornwort orthogroups are all genomes produced in this project, plus A. agrestis 'Bonn' (Li et al. 2020) and A. angustus (Zhang et al. 2020).

Bryophyte orthogroups are all in hornwort analysis plus:

| Genome | Source | Publication |
|--------|--------|-------------|
| Ceratodon purpureus GG1 v1.1 | Phytozome v13 | Carey et al. 2021 |
| Ceratodon purpureus R40 v1.1 | Phytozome v13 | Carey et al. 2021 |
| Entodon seductrix | https://doi.org/10.6084/m9.figshare.17097077 | Yu et al. 2022 |
| Fontinalis antipyretica | http://gigadb.org/dataset/100748 | Yu et al. 2020 |
| Hypnum curvifolium | https://doi.org/10.6084/m9.figshare.17097077 | Yu et al. 2022 |
| Lunularia cruciata 1F | CoGe <sup>1</sup> | Linde et al. 2023 |
| Lunularia cruciata 1M | CoGe <sup>1</sup> | Linde et al. 2023 |
| Marchantia polymorpha v3.1 | Phytozome v13 | Bowman et al. 2017 |
| Physcomitrium patens v3.3 | Phytozome v13 | Lang et al. 2018 |
| Sphagnum fallax v1.1 | Phytozome v13 | Healey et al. 2023 |
| Syntrichia caninervis v2 | CoGe | Silva et al. 2021 |

<sup>2</sup> Translated CDS FASTA was downloaded from CoGe. Gene IDs in sequence headers suggest alternative transcripts are not included. BUSCO results also support that these are primary transcripts so no further processing was done.
<sup>1</sup> Translated CDS FASTA downloaded from CoGe appears to contain alternative transcripts based on sequence headers and BUSCO results. Some transcripts are also translated in multiple frames. The transcripts labelled ".1" and "frame0" were selected for further analysis. The total complete BUSCOs is identical between source and filtered data, but with 15% vs 1% duplicated, respectively.

<details>
<summary>Example procedure for filtering L. cruciata:</summary>

```shell
for i in {00000..23502} ; do grep $i"\.1|" Lunularia_cruciata_1M_64631-CDS-prot.fasta & done | grep "frame0" | getFromFasta.py Lunularia_cruciata_1M_64631-CDS-prot.fasta - > Lunularia_cruciata_1M_64631-CDS-prot.fasta.1transcripts_frame0.fasta
```
</details>
<details>
<summary>All CoGe files had to be renamed using the unique gene identifier:</summary>

```shell
grep ">" Syntrichia_caninervis_65083-CDS-prot.fasta | awk -F"|" '{print $0"\t"$9}' > Syntrichia_caninervis_65083-CDS-prot.fasta.rename-table.tsv
sed -i 's/>//' Syntrichia_caninervis_65083-CDS-prot.fasta.rename-table.tsv
renameFasta.py Syntrichia_caninervis_65083-CDS-prot.fasta Syntrichia_caninervis_65083-CDS-prot.fasta.rename-table.tsv
```
</details>
<br>
Viridiplantae orthogroups include all in bryophyte analysis plus:

| Genome | Source | Publication |
|--------|--------|-------------|
| **Chlorophyte Algae**  |        |
| Chara braunii <sup>*</sup> | OrcAE (https://bioinformatics.psb.ugent.be/orcae/overview/Chbra) | Nishiyama et al. 2018 |
| Chylamydomonas reinhardtii 281 v5.6 | Phytozome v13 | Merchent et al. 2007 |
| Dunaliella salina v1.0 | Phytozome v13 | Polle et al. 2020 |
| Klebsormidium nitens v1.1| http://www.plantmorphogenesis.bio.titech.ac.jp/~algae_genome_project/klebsormidium/   | Hori et al. 2014 |
| Mesotaenium endlicherianum | https://doi.org/10.6084/m9.figshare.9911876.v1 | Cheng et al. 2019 |
| Ostreococcus lucimarinus v2.0 | Phytozome v13 | Palenik et al. 2007 |
| Spirogloea muscicola<sup>^</sup> | https://doi.org/10.6084/m9.figshare.9911876.v1 | Cheng et al. 2019 |
| Volvox carteri v2.1 | Phytozome v13 | Prochnik et al. 2010 |
| Zygnema circumcarinatum SAG 698-1b | https://bcb.unl.edu/Zygnema_4genomes/SAG698-1b/ | Feng et al. 2023 |
| **Lycophytes** |    |
| Diphasiastrum complanatum v3.1 | Phytozome (restricted - in review) |
| Isoetes taiwanensis v1 <sup>*</sup>| CoGe | Wickell et al. 2021 |
| Selaginella moellendorffii v1.0 | Phytozome v13 | Banks et al. 2011 |
| **Pteridophytes**  |        |
| Adiantum capillus-veneris<sup>+</sup> | FernBase (fernbase.org) | Fang et al. 2022 |
| Alsophila spinulosa v3.1 <sup>+</sup> | FernBase (fernbase.org) | Huang et al. 2022 |
| Azolla filiculoides v1.1 | FernBase (fernbase.org) | Li et al. 2018 |
| Ceratopteris richardii v2.1 | Phytozome v13 | Marchant et al. 2022 |
| Marsilea vestita v3 | FernBase (fernbase.org)  | Rahmatpour et al. 2023 |
| Salvinia cucullata v1.2 | FernBase (fernbase.org) | Li et al. 2018 |
| **Gymnosperms** |    |
| Cycas panzhihuaensis | https://db.cngb.org/codeplot/datasets/PwRftGHfPs5qG3gE | Liu et al. 2022 |
| Ginkgo biloba v2 | Ginkgo DB (ginkgo.zju.edu.cn/genome) | Liu et al. 2021 |
| Gnetum montanum | https://datadryad.org/stash/dataset/doi:10.5061/dryad.ht76hdrdr | Wan et al. 2021 |
| Sequoia sempervirens v2.1 | TreeGenesDB | Scott et al. 2020 |
| Welwitschia mirabilis | https://datadryad.org/stash/dataset/doi:10.5061/dryad.ht76hdrdr | Wan et al. 2021 |
| **Angiosperms** |
| Amborella trichopoda v1.0 | Phytozome v13 | Albert et al. 2013 |
| Arabidopsis thaliana TAIR10 | Phytozome v13 | Lamesch et al. 2011 |
| Cinnamomum kanehirae v3 | Phytozome v13 | Chaw et al. 2019 |
| Medicago truncatula Mt4.0v1 | Phytozome v13 | Tang et al. 2014 |
| Nymphaea colorata v1.2 | Phytozome v13 | Zhang et al. 2019 |
| Oryza sativa v7.0 | Phytozome v13 | Ouyang et al. 2007 |
| Solanum lycopersicon ITAG4.0 | Phytozome v13 | Hosmani et al. 2019 |
| Zea mays RefGen v4 | Phytozome v13 | NA |

<sup>*</sup> Translated CDS FASTA appears to contain alternative transcripts, which were removed before analysis
<sup>+</sup> 415 sequences in Alsophila spinulosa and 361 sequences in Adiantum capillus-venerus contained in-frame stop codons. These sites were converted to 'X' in order to include those genes in the analysis
<sup>^</sup> BUSCO results showed a very high percentage of duplicated genes, but analysis of spatial coordinates in the FASTA sequence metadata shows this is not due to the presence of alternative transcripts

<details>
<summary>Filtering Isoetes taiwanensis</summary>

```shell
grep ">" Isoetes_taiwanensis_61511-CDS-prot.fasta | awk -F"|" '{print $0"\t"$9}' > Isoetes_taiwanensis_61511-CDS-prot.fasta.rename-table.tsv
sed -i "s/>//" Isoetes_taiwanensis_61511-CDS-prot.fasta.rename-table.tsv
renameFasta.py Isoetes_taiwanensis_61511-CDS-prot.fasta Isoetes_taiwanensis_61511-CDS-prot.fasta.rename-table.tsv
grep "\-RA" Isoetes_taiwanensis_61511-CDS-prot_renamed.fasta | getFromFasta.py Isoetes_taiwanensis_61511-CDS-prot_renamed.fasta - > Isoetes_taiwanensis_61511-CDS-prot_renamed.primary_transcripts.fasta
```

</details>

<details>
<summary>Filtering Chara braunii</summary>

```shell
removeAlternativeTranscripts.py chbra_iso_noTE_23546_pep.fasta
```
</details>

<details>
<summary>Converting in-frame stop codons (Python)</summary>

```python
import re
with open("Alsophila_spinulosa_v3.1_protein.fa", "r") as infile, open("Alsophila_spinulosa_v3.1_protein.IFSconverted.fa", "w") as outfile:
	for line in infile:
		if line.startswith(">"):
			outfile.write(line)
		else:
			outfile.write(re.sub("\.", "X", line))
```
</details>

<details>
<summary>Check Spirogloea muscicola for alternative transcripts</summary>

Convert FASTA header information to BED format:
```shell
grep ">" Spirogloea_muscicola_gene.pep.fasta | awk -F" |:|=" '{print $6"\t"$7"\t"$8"\t"$1"\t0\t"$9}' | sort -k 1,1 -k 2,2n > Spirogloea_muscicola_gene.pep.fasta.sorted.bed
```
Use bedtools to merge overlapping or touching regions:
```shell
bedtools merge -i Spirogloea_muscicola_gene.pep.fasta.sorted.bed -d 0 -s -c 1,4 -o count,distinct | awk -F"\t" '{ if ( $4 > 1 ) print $0}'
```
</details>

### Orthogroup Inference
```shell
# OrthoFinder version 2.5.4
conda activate orthofinder
orthofinder -M msa -t 24 -S diamond -A mafft -T fasttree -X -f /mnt/data/ps997/projects/hornwort_genomes/analysis/orthofinder/hornworts-only/

```

#### Orthogroups on Hornwort Accessory/Sex Chromosomes
Orthogroups were filtered to find those containing genes on an accessory/sex chromosome (AnagrOXF.S6, Ledus.S6, Noorb.S5, Papea.S5, Papea.S6, PhPhy.S5) in more than one species.

<details>
<summary>Filtering orthogroups and plotting</summary>

```shell
grep "AnagrOXF.S6 Orthogroups.tsv  > Orthogroups_any_acc-sex_chromosomes.AnagrOXF.tsv
grep "Ledus.S6" Orthogroups.tsv  > Orthogroups_any_acc-sex_chromosomes.Ledus.tsv
grep "Noorb.S5" Orthogroups.tsv  > Orthogroups_any_acc-sex_chromosomes.Noorb.tsv
grep "Papea.S5\|Papea.S6" Orthogroups.tsv  > Orthogroups_any_acc-sex_chromosomes.Papea.tsv
grep "Phphy.S5" Orthogroups.tsv > Orthogroups_any_acc-sex_chromosomes.Phphy.tsv

# Create list of orthogroups from each file to count number of occurrences and select only those in more than one species.
cut -f 1 Orthogroups_any_acc-sex_chromosomes.*.tsv > Orthogroups_any_acc-sex_chromosomes.allOGs.txt
pyTable.py Orthogroups_any_acc-sex_chromosomes.allOGs.txt | awk -F"\t" '{ if ($2 > 1) print $1}' > Orthogroups_any_sex_acc_chromosome.non-species-specific.txt

# Manually record presence/absence data for each OG in the species-specfic files in a CSV file
cat Orthogroups_any_sex_acc_chromosome.non-species-specific.txt | while read line ; do grep $line Orthogroups_any_sex_acc_chromosome.*.tsv | cut -f 1 ; done
```

Manually recorded data in CSV format: <a href ="files/orthogroups_acc-sex_chromosomes.csv">orthogroups_acc-sex_chromosomes.csv</a>

```
library(UpSetR)
ogs <- read.csv("orthogroups_acc-sex_chromosomes.csv", header = T)
upset(ogs, sets = c("Phphy", "Papea", "Noorb", "Ledus", "AnagrOXF"), order.by = c("freq", "degree") , mainbar.y.label = "Orthogroups", sets.x.label = "Orthogroups", keep.order =  TRUE)
```

</details>

### Genespace
Primary transcript protein FASTAs were downloaded from HornwortBase/Phytozome. For Phytozome FASTAs, sequences had to be renamed to match those in the GFF. GFF annotation files were downloaded from HornwortBase/Phytozome, then converted into a simplified 4-column BED format. FASTA and BED files were renamed to match each other.

<details>
<summary>Convert Phytozome FASTAs</summary>

```
grep ">" CpurpureusGG1_539_v1.1.protein_primaryTranscriptOnly.fa | awk -F" |=" '{print $0"\t"$7}' > Ceratodon_purpureus_GG1.convert.tsv
renameFasta.py CpurpureusGG1_539_v1.1.protein_primaryTranscriptOnly.fa Ceratodon_purpureus_GG1.convert.tsv
mv CpurpureusGG1_539_v1.1.protein_primaryTranscriptOnly_renamed.fa Ceratodon_purpureus_GG1.fa
```
</details>

<details>
<summary>Convert GFF to simplified BED</summary>

For Phytozome:
```shell
grep -w "gene" Ppatens_318_v3.3.gene.gff3 | gff2bed | awk -F"\t|=|;" '{print $1"\t"$2"\t"$3"\t"$13}' > Ppatens_318_v3.3.gene.bed  
```
For HornwortBase:
```shell
grep -w "gene" Anthoceros_agrestis_Oxford_gene_annotations.gff | gff2bed | awk -F"\t|=|;" '{print $1"\t"$2"\t"$3"\t"$11}' > Anthoceros_agrestis_Oxford.bed
```
</details>

<details>
<summary>Run Genespace</summary>

```shell
conda activate orthofinder
```
In R:
```
library(GENESPACE)
gpar <- init_genespace(wd = "/media/peter/WD14TB-2/projects/hornwort_genomes/analysis/synteny/genespace", path2mcscanx = "~/bin/MCScanX/"
gpar <- run_genespace(gsParam = gpar)
```
</details>


### GO Enrichment of Accessory/Sex Chromosomes
GO enrichment done using goatools, filtered to Bonferroni-corrected p-value <0.05, and plotted in R.

File formatting:
```
# Background genes:
awk -F"\t" '{print $1}' Anthoceros_agrestis_Oxford_gene_functional_annotations.tsv > AnagrOXF_background.txt

# Target genes:
grep "AnagrOXF.S6" AnagrOXF_background.txt > AnagrOXF.S6.txt

# Association file (mapping gene ID to GO terms):
awk -F"\t" '{print $1"\t"$10}' Anthoceros_agrestis_Oxford_gene_functional_annotations.tsv > AnagrOXF_association.txt
sed -i 's/,/;/g' AnagrOXF_association.txt
```

goatools:
```
find_enrichment.py --obo /mnt/data/ps997/resources/GO_databases/go-basic.obo --annofmt id2gos  --pval 0.05 --outfile AnagrOXF.S6.GOenrichment.xlsx AnagrOXF.S6.txt AnagrOXF_background.txt AnagrOXF_associations.txt
```

R filtering and plotting:
```
library(ggplot2)
anagroxfS6 <- read.csv("AnagrOXF.S6.GOenrichment.csv", header = T)
anagroxfS6Filt <- anagroxfS6[anagroxfS6$p_bonferroni < 0.05, ]
anagroxfS6Plot <- ggplot(data = anagroxfS6Filt, aes(x=reorder(anagroxfS6Filt$name, anagroxfS6Filt$study_count), anagroxfS6Filt$study_count, fill = anagroxfS6Filt$NS)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw() +
  xlab("GO Term Description") +
  ylab("Count") +
  labs(fill = "GO Categories") +
  scale_y_continuous(breaks = seq(0, max(anagroxfS6Filt$study_count), 10))
```

<br>
<br>
<br>
<br>



Orphan genes on sex/microchromosomes
RNA replicates for Anthoceros and Notothylas
r8s time tree - R package?
