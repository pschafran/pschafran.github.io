---
layout: post
title:  "Pan-phylum genomes of hornworts revealed conserved autosomes but dynamic accessory and sex chromosomes"
date:   2024-04-30 15:37:17 -0400
categories: jekyll update
---
<section>
<p>Intial Post Date: 2024-04-30</p>
</section>

<section id="main">
<h2>Table of Contents</h2>

1. <a href = "#abstract">Abstract</a>
2. <a href = "#methods">Methods</a>
3. <a href = "#results-and-discussion">Results and Discussion</a>
4. <a href = "#supplementary-info">Supplementary Info</a>

<section id="abstract">
<h2>Abstract</h2>
</section>

<section id="methods">
<h2>Methods</h2>


<section id="assembly">
<h3>1. Assembly</h3>
<p>ONT reads less than 5 kbp were removed and the remainder were assembled with Flye v2.9:</p>


<code>
awk 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 5000) {print header, seq, qheader, qseq}}' < your.fastq > filtered.fastq

flye --nano-hq filtered.fastq -t 24 -o flye 
</code>


Contigs were corrected with Illumina DNA sequence data using Pilon v1.24 in three iterations, with the Pilon output as input each successive round:</p>


<code>
bwa index assembly.fasta
bwa mem -t 24 assembly.fasta Illumina_reads_R1.fq Illumina_reads_R2.fq | samtools sort -o illumina.bam
minimap2 -t 24 assembly.fasta ONT_reads.fq | samtools sort -o ont.bam
java -Xmx50G -jar pilon-1.24.jar --genome assembly.fasta --frags illumina.bam --nanopore ont.bam --output pilon
</code>


<p>Three rounds of Pilon generally corrected >90% of all changes with diminishing returns (and possible over-polishing errors) with further rounds.</p>
<p>HiC libraries were prepared, sequenced, and scaffolded by Phase Genomics (Seattle, WA). TGS-Gapcloser was used to fill gaps between scaffolds with ONT reads and polish filled gaps with Illumina data: </p>

<code>
convertFastqToFasta.py ONT_reads.fq
cat Illumina_reads_R1.fq Illumina_reads_R2.fq > Illumina_reads_combined.fq
tgsgapcloser --scaff scaffolded_assembly.fasta \
    --reads ONT_reads.fasta \
    --ouput assembly.gapclosed
    --ngs Illumina_reads_combined.fq \
    --pilon /home/ps997/bin/pilon-1.24.jar \ # Change to match your system
    --samtools /usr/local/bin/samtools \ # Change to match your system
    --java /usr/bin/java # Change to match your system
</code>

</section> <!--Assembly end-->

<section id="repeat-annotation">
<h3>2. Repeat Annotation</h3>
<p>Repeats were identified in a first pass by EDTA v2:</p>

<code>
EDTA.pl --sensitive 1 --anno 1 --evaluate 1 -t 12 \
--genome scaffolded_assembly.gapclosed.scaff_seqs.fasta \
--repeatmasker /home/ps997/bin/RepeatMasker/RepeatMasker \ # Change to match your system
--cds ~/HornwortBase_20210503/SPECIES_CDS.fna \ # OPTIONAL: I used these because I already had predicted CDS sequences from a previous round of annotation
&> edta.out # Capture the output for potential debugging
</code>

<p>NOTE: EDTA doesn't like long sequence names, you might have to rename them first with some simple Python code:</p>

<code>
infile = open("PGA_assembly.gapclosed.scaff_seqs","r")
outfile = open("PGA_assembly.gapclosed.scaff_seqs.renamed.fasta","w")
counter = 1
for line in infile:
	if line.startswith(">"):
		outfile.write(">%s\n" % counter)
		counter += 1
	else:
		outfile.write("%s" % line)
</code> 

<p>OPTIONAL: EDTA is mostly used as a stand-alone tool these days, but some additional processing may improve the transposable element (TE) library. These steps aim to recover protein-coding genes that were misidentified as LTRs, following steps from the <a href="https://weatherby.genetics.utah.edu/MAKER/wiki/index.php/Repeat_Library_Construction-Advanced">MAKER wiki</a>. Search databases are linked at the bottom of that page.</p>

<code>
# Pull out LTRs not identified as Copia or Gypsy type from the EDTA TE library
extract_unknownLTR.py PGA_assembly.gapcloser.scaff_seqs.renamed.fasta.mod.EDTA.TElib.fa

# BLAST search the unknown LTRs against a curated database of transposons
blastx -query PGA_assembly.gapcloser.scaff_seqs.renamed.fasta.mod.EDTA.TElib.fa.LTRlib.unknown.fa -db Tpases020812 -evalue 1e-10 -num_descriptions 10 -out PGA_assembly.gapcloser.scaff_seqs.renamed.fasta.mod.EDTA.TElib.fa.LTRlib.unknown.fa.blastx -num_threads 12

# Parse the BLAST results and identify unknown LTRs with hits to known transposons
perl Custom-Repeat-Library/transposon_blast_parse.pl --blastx PGA_assembly.gapcloser.scaff_seqs.renamed.fasta.mod.EDTA.TElib.fa.LTRlib.unknown.fa.blastx --modelerunknown PGA_assembly.gapcloser.scaff_seqs.renamed.fasta.mod.EDTA.TElib.fa.LTRlib.unknown.fa

# Combine the LTRs that have BLAST hits to the main LTR library file. Rename the unknown LTRs file. 
cat PGA_assembly.gapcloser.scaff_seqs.renamed.fasta.mod.EDTA.TElib.fa.LTRlib.known.fa identified_elements.txt > PGA_assembly.gapcloser.scaff_seqs.renamed.fasta.mod.EDTA.TElib.fa.LTRlib.known.final.fa
mv unknown_elements.txt PGA_assembly.gapcloser.scaff_seqs.renamed.fasta.mod.EDTA.TElib.fa.LTRlib.unknown.final.fa

# BLAST search the LTR library against UNIPROT plant protein database. 
blastx -db uniprot_sprot_plants.fasta -query PGA_assembly.gapcloser.scaff_seqs.renamed.fasta.mod.EDTA.TElib.fa.LTRlib.known.final.fa -out uniprot_plant_blast.out -num_threads 12

# Parse the BLAST results and remove LTRs that matched plant proteins.
perl ProtExcluder1.1/ProtExcluder.pl uniprot_plant_blast.out PGA_assembly.gapcloser.scaff_seqs.renamed.fasta.mod.EDTA.TElib.fa.LTRlib.known.final.fa

# Mask the repeats in the genome using this new LTR library. The `--xsmall` option softmasks the genome, which is preferred by BRAKER. 
/home/ps997/bin/RepeatMasker/RepeatMasker -noisy -a -gff -u -pa 24 --xsmall -lib PGA_assembly.gapcloser.scaff_seqs.renamed.fasta.mod.EDTA.TElib.fa.LTRlib.known.final.fanoProtFinal PGA_assembly.gapcloser.scaff_seqs.renamed.fasta

# Calculate stats
RepeatMasker/util/buildSummary.pl -useAbsoluteGenomeSize PGA_assembly.gapcloser.scaff_seqs.renamed.fasta.out > PGA_assembly.gapcloser.scaff_seqs.renamed.fasta.repeat-summary.txt
</code>

<h4>2a. Tandem Repeats Finder</h4>
<p>The software `tandem repeats finder` (TRF) can find additional tandem repeats that are associated with centromeres and telomeres.</p>

<p>Run TRF:</p>

<code>
trf Phymatoceros_genome_Phase_scaffolded.fasta 2 7 7 80 10 50 2000 -h -d -m -ngs > trf_out.txt
</code>

<p>Convert TRF output to GFF and BED formats:</p>

<code>
trf2gff.py trf_out.txt 50 > trf_out_min50.gff
awk -F"\t" '{print $1"\t"$2"\ttandem_repeat\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' trf_out_min50.gff > trf_out_min50.renamed.gff
gff2bed < trf_out_min50.gff > trf_out_min50.bed
</code>

<p>TRF will find elements already identified by EDTA. To filter those out, we can remove the overlapping items:</p>

<code>

</code>

</section> <!--Repeat annotation end-->

<section id="gene-prediction">
<h3>3. Gene Prediction</h3>
<p>Gene models were predicted using BRAKER (), with input consisting of Illumina RNA reads mapped to the softmasked genome using HISAT2 () and predicted hornwort proteins from published *Anthoceros* genomes (Li et al 2020, Zhang et al. 2020). BRAKER output files were screened for genes with in-frame stop codons, which were marked as pseudogenes in the corresponding GTF file. Genes were renamed to contain their respective scaffold/contig name plus a number incremented by 100, restarting at the beginning of each scaffold/contig. Subsets of primary transcripts were created by selecting the longest transcript associated with each gene.
</p><br>
<b>Code for mapping RNA reads and running BRAKER and processing output</b><br>

<code>shell
hisat2-build -p 12 PGA_assembly.gapcloser.scaff_seqs.renamed.masked_assembly.fasta PGA_assembly.gapcloser.scaff_seqs.renamed.masked_assembly.fasta
hisat2 -p 12 -x PGA_assembly.gapcloser.scaff_seqs.renamed.masked_assembly.fasta -1 RNA_R1.fq.gz -2 RNA_R2.fq.gz 2> hisat-align.out | samtools sort -o PGA_assembly.gapcloser.scaff_seqs.renamed.masked_assembly.RNAmapped.bam
</code>

<code>
braker.pl \
--genome PGA_assembly.gapcloser.scaff_seqs.renamed.masked_assembly.fasta \
--bam PGA_assembly.gapcloser.scaff_seqs.renamed.masked_assembly.RNAmapped.bam \
--prot_seq Hornwort_orthogroups.faa \
--prg=gth --gth2traingenes \
--verbosity 3 \
--cores 12 \
--nocleanup \
--softmasking \
</code>

<b>Filter genes with in-frame stop codons</b><br>
<p>BRAKER will sometimes predict proteins that contain in-frame (internal) stop codons. In the BRAKER-produced CDS/AA FASTA files, the sequences in the bad regions are masked with N (CDS) or X (AA), but the GTF file will still contain annotation that creates a bad sequence. Following NCBI protocol for genes that are 'broken' but are not thought to be pseudogenes, these will get annotated with `pseudo=true`.
The easiest starting point is the `bad_genes.lst` if you ran BRAKER with the `--nocleanup` option. If you didn't, make a new translation from the GTF:</p>

<code>
gffread -y proteins.fasta -g genome.fasta augustus.hints.gtf
</code>

Then search for sequences with periods (representing stop codons), and write them to `bad_genes.lst`

<code>
with open("proteins.fasta", "r") as infile, open("bad_genes.lst", "w") as outfile:
	bad_genes = []
	for line in infile:
		if line.startswith(">"):
			seqid = line.strip(">|\n")
		elif "." in line and seqid not in bad_genes:
			bad_genes.append(seqid)
			outfile.write("%s\n" % seqid)
</code>

<b>Rename Contigs and Genes</b>
<p>You'll probably want to rename the contigs and genes in the fasta and gff/gtf files associated with each genome. The input genome and any annotations must match (e.g. can't use RepeatModeler annotations if you renamed the genome for EDTA). Here it may be advisable to add a unique ID or version number to link this particular genome assembly and annotation. Comment lines can be added to the header of the final GFF/GTF and fasta metadata can be added to sequence names. In this example, I start with the genome that was temporarily renamed to run EDTA.</p>

<b>Naming conventions for hornwort genomes:</b>
<ol>
	<li>Scaffolds: Two letters of genus + three letters of species (optional information after this e.g cultivar, sex) + period + S/C (for scaffold or contig) + number (order determined by longest to shortest length).</li>
	<li>Genes: Scaffold ID + G + six digit number unique to each gene. Genes increment by 100 and numbering restarts for each sequence.</li>
	<li>Transcripts: Gene ID + "." + transcript number (increments by 1).</li>
</ol>

<code>bash
# To hash the date for a unique ID
md5sum <(date)
6cde96438c2e713efa5c285e4fe3a62d


# Get the length of all sequences in order to rename from shortest to longest
getFastaSeqLengths.py PGA_assembly.gapcloser.scaff_seqs.renamed.fasta
sort -k2,2nr PGA_assembly.gapcloser.scaff_seqs.renamed.fasta_sequence_lengths.tmp > genome.fa_sequence_lengths.tsv.sorted.tmp

# Change "j" to the part of the names that will be the same in all sequences. Change "id" if using a unique ID for this genome and its annotations
awk -F"\t" -v i="1" -v j="AnagrOXF.C" -v id="6cde96438c2e713efa5c285e4fe3a62d" '{ print $1"\t"j""i++" id="id }' genome.fa_sequence_lengths.tsv.sorted.tmp > genome.fa_sequence_lengths.tsv.new_contig_names.tsv

### NOTE: Here I manually edit the *new_contig_names.tsv file to change scaffold/contig designations
renameFastaAndReorder.py PGA_assembly.gapcloser.scaff_seqs.renamed.fasta genome.fa_sequence_lengths.tsv.new_contig_names.tsv

# Rename a BRAKER GTF file
renameGTF_Phytozome.py -i augustus.hints.gtf --contig-table genome.fa_sequence_lengths.tsv.new_contig_names.tsv --assembly-id 6cde96438c2e713efa5c285e4fe3a62d

</code>



<b>Create pseudo-gene annotated GTF</b><br>
<p>Run `renameGTF_Phytozome.py` with `--bad-genes bad_genes.lst`. It will add `pseudo=true` to the GTF file to mark broken genes in [NCBI style](https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/).</p>

</section> <!--Gene prediction end-->

<section id="synteny">
<h3>Synteny</h3>

<b>Extracting syntenic blocks</b>
<p>Information about syntenic blocks was extracted from GENESPACE output file in `GENESPACE-DIR/results/syntenicBlock_coordinates.tsv`.</p>

<p>To get stats about block size between pairs of genomes:</p>

<code>
grep "BONN" syntenicBlock_coordinates.csv | grep "Oxford" | sed 's/,/\t/g' | cut -f 5,12 | sort | uniq | cut -f 2 | summaryStats.R
</code>

<p>Note that the sort and uniq commands are used to avoid double-counting the same block in both orientations.</p>

</section> <!--Synteny end-->

</section> <!--Methods end-->

<section id="results-and-discussion">
<h2>Results and Discussion</h2>
</section> <!--Results-discussion end>



<section id="supplementary-info">
<h2>Supplementary Info</h2>
</section>

</section> <!--Main end-->