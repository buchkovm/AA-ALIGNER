AA-ALIGNER: Allele-Aware ALignments for Investigating GeNetic Effects on Regulation

Enrichment of one allele in quantitative sequence data, or allelic imbalance, can indicated allelic differences in transcription factor binding and/or gene transcriptional regulation.  Reference mapping biases at heterozygous sites hinder accurate identification of allelic imbalance. We have developed the AA-ALIGNER pipeline to overcome reference mapping biases and identify allelic imbalance in quantitative sequence data. AA-ALIGNER can accurately identify imbalance using any available or even no sample genotype data.

AA-ALIGNER 1.0 - AA-ALIgNER for a single machine
usage: perl AA_ALIGNER.pl [options]... [fastq files(s)]
GENERAL OPTIONS:
-configuration_file <configuration file>	configuration file (required)
-out_directory	<output directory>	directory for permanent files (default=./)
-temp_directory <temporary directory>	directory for intermediate files (default=./)
-out_prefix <output prefix>		prefix for output bam file and imbalance files

RUN OPTIONS:
-create_reference 			create a sample-specific reference genome (default=OFF)
-filter_fastq	filters fastq files before alignment to convert base qualities to Sanger format, trim sequences (if necessary), and removes low quality sequence reads (default=OFF)
-align					aligns sequences to specified genomes (default=OFF)
-filter_alignment	filters alignments to remove reads that map multiple sites,fall in blacklist regions, and are possible PCR duplicates (default = OFF)
-find_imbalance	dentifies sites allelic imbalance in  sequence alignment (default=OFF)
-all					runs entire alignment pipeline filter_fastq, align,								filter_alignments, and find_imbalance (default=ON)

AA-ALIGNER CLUSTER 1.0 - AA-ALIGNER for a cluster (currently implemented for LSF only)
usage: perl AA_ALIGNER_cluster.pl [options]... [fastq files(s)]
GENERAL OPTIONS:
-configuration_file <configuration file>			configuration file (required)
-cluster_configuration <cluster configuration file>	configuration file for cluster settings
-out_directory	<output directory>			directory for permanent files (default=./)
-temp_directory <temporary directory>			directory for intermediate files (default=./)
-out_prefix <output prefix>				prefix for output bam file and imbalance files

RUN OPTIONS:
-create_reference 	create a sample-specific reference genome (default=OFF)
-filter_fastq	filters fastq files before alignment to convert base qualities to Sanger format, trim sequences (if necessary), and removes low quality sequence reads (default=OFF)
-align	aligns sequences to specified genomes (default=OFF)
-filter_alignment	filters alignments to remove reads that map multiple sites,fall in blacklist regions, and are possible PCR duplicates (default = OFF)
-find_imbalance	identifies sites allelic imbalance in sequence alignment (default=OFF)
-all	runs entire alignment pipeline filter_fastq, align, filter_alignments, and find_imbalance (default=ON)
Example Usage
Create usage examples are given for AA_ALIGNER but also apply to AA_ALIGNER CLUSTER with the addition of the cluster_configuration option

FASTQ Files used in the example can be found in AA_ALIGNER/example/fastq/
Create a custom reference file



System Dependencies
GENERAL
Perl
https://www.perl.org/

Perl Modules:
Getopt::Long
Math::CDF

Java
https://java.com/en/download/

Python
https://www.python.org/

FASTQ FILTERING
FastQC	
http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

GALAXY SOURCE for FASTQ_groomer
https://wiki.galaxyproject.org/Admin/GetGalaxy

FASTX Toolkit
http://hannonlab.cshl.edu/fastx_toolkit/

SEQUENCE ALIGNMENT
GSNAP/GMAP
http://research-pub.gene.com/gmap/

FILTER ALIGNMENT 
Picard Tools
http://broadinstitute.github.io/picard/

SAMtools
http://samtools.sourceforge.net/

FIND ALLELIC IMBALANCE
SAMtools
http://samtools.sourceforge.net/

GSNAP/GMAP
http://research-pub.gene.com/gmap/

