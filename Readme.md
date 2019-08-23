EnsembleSV
=
A general pipeline for ensemble calling *Structural Variations* (SVs), also referred to as *Novel Adjacencies* (NAs) with multiple methods 
in short-read paired-end Illumina/linked sequencing data as well as in 3rd-generation long ONt/PacBio reads.

Supported SV inference methods:

**short reads**
* *SvABA*
* *Manta*
* *Lumpa*

**linked reads**
* *LongRanger*
* *GROCSVs*
* *NAIBR*

**long (ONT/PacBio) reads**
* **Sniffles**
* **PBSV**

Supported sequencing data is basically described by what type of data the SV inference methods support. 
Current supported input data includes:
* Illumina paired-end short reads
* Illumina/10X Genomics linked paired-end short reads
* Oxford Nanopore Technologies long reads
* Pacific Biosciences long reads

Description
= 

Input data is expected to be in a form of aligned bam files. 
We recommend that alignment is performed to the same reference, but differences in chromosomal names are permitted.
**Note**: do not use alignments to different *versions* of the references, as that breaks coordinate concordance.

Pipeline is implements via *snakemake* workflow manager.

Basic workflow can be described as follows:
1. Method-specific SV inference on all method-compatible alignments in the input
2. Merging of sequencing-technology-specific SV calls into *sensitive-technology-SV-sets*
    * merging is performed with **SURVIVOR**
    * if a method produces more than one SV callset (e.g., separate indel and SV sets) they are merged into method-specific callset with **RCK** utilities 
    * only SVs longer than a `min_len` size threshold are retained
    * for short-read SVs only SVs supported by at least 2 alignment-methods are retained
3. Merging with **SURVIVOR** of *sensitive-technology-SV-sets* into *sensitive-sample-SV-set*
4. Filtration of *sensitive-sample-SV-set* into *specific-sample-SV-set*
     

Installation
=

Usage
=

Contribution
=

Issues
=

Citation
=