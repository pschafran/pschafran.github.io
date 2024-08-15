---
layout: post
title:  "Eremid Genomics and Catawba College -- Illumina Next Generation Sequencing Workshop"
---

## Useful Links

1. [Geneious Prime Free Trial](https://www.geneious.com/free-trial)
2. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
3. [Sequence Data Files (Box Drive)](https://cornell.box.com/s/o5y164l09syhttirifhv3ow1qhsht1en)
4. [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/)
5. [GeSeq](https://chlorobox.mpimp-golm.mpg.de/geseq.html)
6. [BV-BRC Taxonomic Classification Service](https://www.bv-brc.org/app/TaxonomicClassification)
7. [CIPRES Science Gateway](https://www.phylo.org)
8. [Galaxy Webserver](https://usegalaxy.org/)
9. [Miniconda Install](https://docs.anaconda.com/miniconda/#quick-command-line-install)


## Unix/Linux Command Line Plastome Assembly

### Software

* [fastp](https://github.com/OpenGene/fastp)
* [SPAdes](https://github.com/ablab/spades)
* [getOrganelle](https://github.com/Kinggerm/GetOrganelle)
* [NOVOPlasty](https://github.com/ndierckx/NOVOPlasty)
* [iqtree](https://github.com/iqtree/iqtree2)

### Install

Install conda if you haven't already. Following are instructions for Linux and Mac -- pick correct one for your computer. Follow the on-screen prompts using default settings.

```
# MacOS ARM CPU
mkdir -p ~/miniconda3
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh -o ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh

# MacOS Intel CPU
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o ~/miniconda3/miniconda.sh

# Linux
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
```

Intialize

```
~/miniconda3/bin/conda init bash
~/miniconda3/bin/conda init zsh
```

Restart your terminal session for it to take effect.

```
conda install -c bioconda -c conda-forge fastp spades getorganelle novoplasty iqtree matplotlib
```

### Inspect Data

Trim adapters, discard low quality reads with fastp. 

```
fastp -i R1.fastq.gz -I R2.fastq.gz -o R1.fastp.fastq.gz -O R2.fastp.fastq.gz -5 -3 --detect_adapter_for_pe -p --html SAMPLE.html --json SAMPLE.json
```

Inspect fastp html file, check for any abnormalities. 

(Optional) Taxonomic classification to check for bacterial contaminants. https://www.bv-brc.org/app/TaxonomicClassification

### Method 1: Reference-based filtering followed by de novo assembly

Get a complete reference sequence (FASTA format) from NCBI GenBank or other source. 

Index reference sequence with Bowtie2.

```
bowtie2-build reference.fasta reference.fasta
``` 

Map reads to reference, saving mapped reads to new files. Adjust `-p` to match your computer's # CPUs. Using `--very-sensitive-local` enhances short matches for more distantly related species.

```
bowtie2 -p 12 -x reference.fasta -1 R1.fastp.fastq.gz -2 R2.fastp.fastq.gz --al-conc-gz R%.aln.fastq.gz --very-sensitive-local &> /dev/null
```

Assemble mapped reads with SPAdes. Change `-t` to match your computers # of CPUs. Change `-o` to whatever the sample name is. 

```
spades.py -1 R1.aln.fastq.gz -2 R2.aln.fastq.gz -t 12 -o SAMPLE_NAME --careful
```

Check the `scaffolds.fasta` file in the output folder. If successful, there should be three large scaffolds, approx. 100k, 25k, and 15k in length, representing the LSC, SSC, and IR regions. Manually stitch them together in Geneious based on the reference. 
If more and/or smaller scaffolds are produced, try mapping them to the reference in Geneious to see what's missing. 

### Method 2: Completely de novo with getOrganelle

Input all adapter-trimmed reads into getOrganelle. To set up getOrganelle for the first time:

```
# This downloads all databases for different organisms and organelles. See help menu/documentation to download certain ones.
get_organelle_config.py -a all
```
Run getOrganelle for an embryophyte plant plastid. See documention for other organism and organelle options.

```
get_organelle_from_reads.py -1 R1.fastp.fastq.gz -2 R2.fastp.fastq.gz -F embplant_pt -t 12
```

Inspect output files. Hopefully you'll get a complete circular sequence with filename like: `embplant_pt.K115.complete.graph1.1.path_sequence.fasta`. 

Remap reads to assembled plastome to evaluate. 

```
evaluate_assembly_using_mapping.py -f embplant_pt.K115.complete.graph1.1.path_sequence.fasta -1 R1.fastp.fastq.gz -2 R2.fastp.fastq.gz -t 12 -o evaluate_output --draw
```

Check coverage levels plotted in `mapping.pdf`. "Matched"  level should be pretty even (though noisy) across LSC/SSC, maybe be different in IRs.

### Method 3: Completely de novo with NOVOPlasty.

NOVOPlasty requires a starting sequence (FASTA format) that closely matches your organism. Find a sequence for a conserved gene (e.g. rbcL) from GenBank or another source. Copy the file into your working directory. 

//TODO
