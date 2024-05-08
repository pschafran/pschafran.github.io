---
layout: post
title:  "Pan-phylum genomes of hornworts revealed conserved autosomes but dynamic accessory and sex chromosomes"
date:   2024-04-30 15:37:17 -0400
categories: jekyll update
---

<h2>Table of Contents</h2>
1. <a href = "#abstract">Abstract</a>
2. <a href = "#introduction">Introduction</a>
3. <a href = "#results-and-discussion">Results and Discussion</a>
4. <a href = "#methods">Methods</a>
5. <a href = "#supplementary-info">Supplementary Info</a>

<section id="abstract">
<h2>Abstract</h2>
</section>


<section id="introduction">
<h2>Introduction</h2>
</section>


<section id="results-and-discussion">
<h2>Results and Discussion</h2>
</section>


<section id="methods">
<h2>Methods</h2>
</section>


Code for mapping RNA reads and running BRAKER and processing output</summary>

```shell
hisat2-build -p 12 PGA_assembly.gapcloser.scaff_seqs.renamed.masked_assembly.fasta PGA_assembly.gapcloser.scaff_seqs.renamed.masked_assembly.fasta
hisat2 -p 12 -x PGA_assembly.gapcloser.scaff_seqs.renamed.masked_assembly.fasta -1 RNA_R1.fq.gz -2 RNA_R2.fq.gz 2> hisat-align.out | samtools sort -o PGA_assembly.gapcloser.scaff_seqs.renamed.masked_assembly.RNAmapped.bam
```
```
braker.pl \
--genome PGA_assembly.gapcloser.scaff_seqs.renamed.masked_assembly.fasta \
--bam PGA_assembly.gapcloser.scaff_seqs.renamed.masked_assembly.RNAmapped.bam \
--prot_seq Hornwort_orthogroups.faa \
--prg=gth --gth2traingenes \
--verbosity 3 \
--cores 12 \
--nocleanup \
--softmasking \
```
**Filter genes with in-frame stop codons**<br>
BRAKER will sometimes predict proteins that contain in-frame (internal) stop codons. In the BRAKER-produced CDS/AA FASTA files, the sequences in the bad regions are masked with N (CDS) or X (AA), but the GTF file will still contain annotation that creates a bad sequence. Following NCBI protocol for genes that are 'broken' but are not thought to be pseudogenes, these will get annotated with `pseudo=true`.
The easiest starting point is the `bad_genes.lst` if you ran BRAKER with the `--nocleanup` option. If you didn't, make a new translation from the GTF:

```
gffread -y proteins.fasta -g genome.fasta augustus.hints.gtf
```

Then search for sequences with periods (representing stop codons), and write them to `bad_genes.lst`

```
with open("proteins.fasta", "r") as infile, open("bad_genes.lst", "w") as outfile:
	bad_genes = []
	for line in infile:
		if line.startswith(">"):
			seqid = line.strip(">|\n")
		elif "." in line and seqid not in bad_genes:
			bad_genes.append(seqid)
			outfile.write("%s\n" % seqid)
```

<section id="supplementary-info">
<h2>Supplementary Info</h2>
</section>
