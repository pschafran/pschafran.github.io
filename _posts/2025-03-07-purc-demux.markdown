---
layout: post
title:  "Oxford Nanopore Amplicon Sequencing and OTU Inference with PURC"
date:   2025-03-07 10:00:00 -0400
---
<section>
<p>Last Updated: 2025-03-07</p>
<p><a href="https://github.com/pschafran/purc/tree/demux-input">Github Repo</a></p>
</section>

<section>
<h2>Intro</h2>
<p>Sequencing of complex amplicon populations is widely used in metabarcoding and polyploid systematics. While Illumina and PacBio HiFi are the preferred data types due to their higher accuracy,
Oxford Nanopore (ONT) fills a useful niche for experiments requiring long read lengths and low throughput. To process ONT data and construct OTUs for further analysis, I modified the PURC program
to handle the ONT data. Below I describe the steps for preprocessing ONT raw data and using the new development version of PURC.
</p>

<h2>General Notes</h2>
<ul>
        <li>IMO quality is more important than quantity in reads, so I employ very stringent filtering: Qscore > 20, barcodes present on both read ends, duplex reads only. 
	If too few reads are produced at the end, these filters can be relaxed, with the risk that less accurate reads can introduce more noise in the final data that complicates interpretation 
	(i.e. too many OTUs).</li>
</ul>

<h2>Dependencies</h2>
<p>See the <a href="https://github.com/pschafran/purc">PURC github</a> for description of install PURC and its dependencies.</p>
<ul>
	<li><a href="https://github.com/nanoporetech/dorado/">dorado</a> is required for basecalling and demultiplexing</li>
	<li><a href="https://github.com/samtools/samtools">samtools</a> is used for converting SAM/BAM files</li>
</ul>

<h2>Working Environment</h2>
<p>I do basecalling and demultiplexing in the sequencing directory immediately above the POD5 files. Demultiplexed FASTQ files get copied to a new working directory for running PURC.   
</p>

<h2>Basecalling</h2>
<p>Select the highest quality basecalling model for the type of data you have. You may need to download the model (`dorado download ...`) if it's your first time running.
</p>

```
dorado duplex --min-qscore 20 dna_r10.4.1_e8.2_400bps_sup\@v5.0.0 pod5_skip/ > duplex.bam
```

<h2>Demultiplexing</h2>
<p>First, select just the duplex reads from the basecalling output (which contains simplex and duplex reads). Then use dorado to detect barcodes and separate sequences into separate files in 
a new directory called `demultiplex`. Select the kit you used for making your library.
</p>

```
samtools view -h --tag dx:1 duplex.bam > duplex_reads.sam

dorado demux -o demultiplex --kit-name SQK-NBD114-24 --emit-fastq --barcode-both-ends duplex_reads.sam
```
<h2>PURC OTU Clustering</h2>
<p>I copy demultiplexed FASTQ files to a new directory, but the analysis could be done in the same location. You can also rename the files at this point -- PURC will treat the filenames as sample 
names in its output. Either FASTQ or FASTA files are accepted (if FASTQ, PURC will convert them to FASTA before processing). 
</p>

```
mkdir -p /data/projects/purc/demultiplex/
cp ./demultiplex/*fastq /data/projects/purc/demultiplex/
```

<p>PURC requires a configuration file that sets all the run parameters. <a href="https://github.com/pschafran/purc/blob/master/purc_configuration_example.txt">See example</a>. 
A new line needs to be added to the config file (anywhere should be fine). `Input_sequence_dir` is the location of the FASTQ files created by `dorado demux`. **All files in this directory will be 
treated as input**. 
</p>

```
Input_sequence_dir = /data/projects/purc/demultiplex/ # Directory containing fasta/fastq files with sequences that are already demultiplexed with barcodes removed
```

<p>A new option for `Mode` is added here to tell PURC to you are using demultiplexed input files.
</p>

```
Mode	= 2
```

<p>Right now, further splitting of each sample into separate loci is not supported, so only a single locus should be listed (same on the primer lines).
</p>

```
[Loci]
Locus_name                              = LFY   # The order must be identical to that in Locus-barcode-taxon_map
Locus-barcode-taxon_map =

[Primers]
Forward_primer                  = GATCTTTATGAACAATGTGGGA        # 5'-3' Order must match locus_name and locus-barcode-taxon_map!
Reverse_primer                  = GAAATACCTGATTTGTAACCTA                # 5'-3' Order must match locus_name and locus-barcode-taxon_map!
```

<p>Only OTU clustering has been tested at this time. Because DADA2 uses Illumina and PacBio-specific error models for ASV inference, and ONT tends to produce different types of errors, testing is
needed to see if DADA2 can be made to work with ONT. 
</p>

<p>The final config file should look like this:
</p>

```
***Important ones***
[Files]
Input_sequence_file             =       # Fasta/fastq file containing your sequencing reads
in_Barcode_seq_file             =               # Fasta file containing barcode names and sequences; the input for making blast database
in_RefSeq_seq_file              = LFY_refs_diploids.fasta       # Fasta file containing reference names and sequences; the input for making blast database
Input_sequence_dir              = ./demultiplexed_files/ # Directory containing fasta/fastq files with sequences that are already demultiplexed with barcodes removed
Output_prefix                   = purc_run                              # The prefix for output files
Output_folder                   = purc_out                              # The folder in where all purc output will be saved
Log_file                                = log                                   # If left empty, default time-stamped log file name will be used

[Loci]
Locus_name                              = LFY   # The order must be identical to that in Locus-barcode-taxon_map
Locus-barcode-taxon_map =

[Primers]
Forward_primer                  = GATCTTTATGAACAATGTGGGA        # 5'-3' Order must match locus_name and locus-barcode-taxon_map!
Reverse_primer                  = GAAATACCTGATTTGTAACCTA                # 5'-3' Order must match locus_name and locus-barcode-taxon_map!

***Setting up PURC run***
[PPP_Configuration]
Mode    = 2             # 0: Check concatemers and then full run
                                                                # 1: Skip concatemer-checking
                                                                # 2: Use files that are already demultiplexed
Multiplex_per_barcode   = 0             # 0: Each barcode contains only one sample
                                                                # 1: Each barcode contains multiple samples
Dual_barcode                    = 0             # 0: Barcodes only on one primer
                                                                # 1: Unique barcodes on both primers
                                                                # 2: Same barcode on both primers
Barcode_detection               = 1             # 0: Search barcode in entire sequences; will produce "ErrMidBC" flag if barcodes are not at the ends of the sequence
                                                                # 1: Search barcode only at the ends of sequences
Recycle_chimeric_seq    = 0             # 0: Do not recycle
                                                                # 1: Split chimeric sequences into respective locus
Recycle_no_barcoded_seq = 0             # 0: Do not recycle
                                                                # 1: Use Smith-Waterman local alignment (more sensitive) to find barcodes in those sequences that BLAST failed to detect barcode
Clustering_method = 1                   # 0: Use DADA2 ASV inference
                                                                # 1: Use Vsearch OTU clustering
                                                                # 2: Use both clustering methods
Align                                   = 1             # 0: No aligning attempted.
                                                                # 1: Final consensus sequences will be aligned with MAFFT

***Miscellaneous***
[DADA Filtering Parameters]
minLen                          = 0     # The minimum length to keep a read (leave as 0 to use Tukey's outlier detection method)
maxLen                          = 0     # The maximum length to keep a read (leave as 0 to use Tukey's outlier detection method)
maxEE                           = 5 #  After truncation, reads with higher than maxEE "expected errors" will be discarded.
Use_OTU_priors = FALSE # Set to TRUE to use OTU output sequences as priors for ASV inference; only compatible with Clustering_method = 2 (default FALSE)

[Clustering Parameters]
clustID1                                = 0.997         # The similarity criterion for the initial VSEARCH clustering
clustID2                                = 0.995         # The similarity criterion for the second clustering
clustID3                                = 0.990         # The similarity criterion for the third clustering
clustID4                                = 0.997         # The similarity criterion for the FINAL clustering to remove identical/near identical clusters
sizeThreshold1                  = 1             # The min. number of sequences/cluster necessary for that cluster to be retained (set to 2 to remove singletons, 3 to remove singletons and doubles, etc)
sizeThreshold2                  = 4                     # The min. number of sequences/cluster necessary for that cluster to be retained (set to 2 to remove singletons, 3 to remove singletons and doubles, etc)

[Chimera-killing Parameters]
abundance_skew                  = 1.9

[Lima Override]
Lima_override = 1 # Set to 1 to use BLAST-based demultiplexing instead of Lima (Linux only)

[Dependencies] # by default assumed to be available in the shell's PATH
Vsearch                                 = vsearch                                                               # The name of (or path to) Vsearch executable
Cutadapt                                = cutadapt                                                              # The name of (or path to) Cutadapt executable
MAFFT                                   = mafft                                                         # The name of (or path to) MAFFT executable
Rscript                                 = Rscript                                                               # The name of (or path to) Rscript executable
blastn                                  = blastn                                                                # The name of (or path to) blastn executable


[Miscellaneous]
Threads                                 = 8                                             # The number of threads for running BLAST
Verbose_level                   = 2                                             # 0: quiet, 1: noticeable (only uchime+cutadapt output), 2: annoying (all usearch+cutadapt+muscle output)
Remove_intermediates    = 0                                             # 0: keep all the intermediate files, 1: remove intermediate files
Barcode_blastDB                 = barcode_blastdb               # The name of the barcode blast database
RefSeq_blastDB                  = refseq_blastdb                # The name of the reference sequence blast database
seq_name_toErase                = m131213_174801_42153_c100618932550000001823119607181400_                      # Remove this from sequence names
```

<h2>Run PURC</h2>
<p>Once the config file is set, run `purc.py` (from the demux-input branch on Github).
</p>

```
purc.py config.txt
```

<p>Information printed to the screen will look like this:
</p>

```
Start Time: 07/03/2025 14:34:58
Output folder: /Users/peter/purc/sandbox/purc_out
Checking dependencies...
Checking demultiplexed input files...
	Found files:
	033e00f178fd3a2bc8d9f3d79d69d63e2faf143d_SQK-NBD114-24_barcode08.fasta
	033e00f178fd3a2bc8d9f3d79d69d63e2faf143d_SQK-NBD114-24_barcode09.fasta
	033e00f178fd3a2bc8d9f3d79d69d63e2faf143d_SQK-NBD114-24_barcode03.fasta
	033e00f178fd3a2bc8d9f3d79d69d63e2faf143d_SQK-NBD114-24_barcode01.fasta
	033e00f178fd3a2bc8d9f3d79d69d63e2faf143d_SQK-NBD114-24_barcode05.fasta
	033e00f178fd3a2bc8d9f3d79d69d63e2faf143d_SQK-NBD114-24_barcode07.fasta
	033e00f178fd3a2bc8d9f3d79d69d63e2faf143d_SQK-NBD114-24_barcode02.fasta
	033e00f178fd3a2bc8d9f3d79d69d63e2faf143d_SQK-NBD114-24_barcode06.fasta
	033e00f178fd3a2bc8d9f3d79d69d63e2faf143d_SQK-NBD114-24_barcode04.fasta
Renaming sequences...
Reorienting sequences based on references...
	Reorienting 225 sequences...
	Writing reoriented sequences to /Users/peter/purc/sandbox/purc_out/tmp/purc_run_1_bc_trimmed.reoriented.fa
Removing primers...
	745 sequences survived after primer-trimming
Annotating seqs...
	745 sequences annotated
	0 sequences cannot be classified
Splitting sequences into a folder/file for each locus...
Clustering/dechimera-izing seqs...

Working on: LFY...


Putting all the sequences together...

OTU Runtime: 2.44 minutes
```

<p>Some downstream steps are not working yet, but you should be able to get the final OTUs for each sample (in `./purc_out/LOCUS/TAXON/` directories) or those sequences combined 
(`./purc_out/LOCUS/LOCUS.fa`).
</p>

</section>
