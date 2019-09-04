EnsembleSV
=
A general pipeline for ensemble calling *Structural Variations* (SVs), also referred to as *Novel Adjacencies* (NAs) with multiple methods 
in short-read paired-end Illumina/linked sequencing data as well as in 3rd-generation long ONt/PacBio reads.

Supported SV inference methods:

**short reads**
* *SvABA*
* *Manta*
* *Lumpy*

**linked reads**
* *LongRanger*
* *GROCSVs*
* *NAIBR*

**long (ONT/PacBio) reads**
* *Sniffles*
* *PBSV*

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
    * merging is performed with *SURVIVOR*
    * if a method produces more than one SV callset (e.g., separate indel and SV sets) they are merged into method-specific callset with *RCK* utilities 
    * only SVs longer than a `min_len` size threshold are retained
    * for short-read SVs only SVs supported by at least 2 alignment-methods are retained
3. Merging with *SURVIVOR* of *sensitive-technology-SV-sets* into *sensitive-sample-SV-set*
4. Filtration of *sensitive-sample-SV-set* into *specific-sample-SV-set*
     

Installation
=

Clone this repository into the **assests** location (i.e., place, where the **reference** version of the workflow resides):

````
git clone https://github.com/aganezov/EnsembleSV.git
```` 

Usage
=
*EnsembleSV* is designed to be utilized one-sample at a time.
We assume that `path` is the location, where the analysis will take place. 
We also assume that `soft/EnsembleSV` is the path for the cloned repository.
Then:
````bash
cd path
mkdir sv && cd sv
ln -s soft/EnsembleSV/*.snakefile .
ln -s soft/EnsembleSV/conda
ln -s soft/EnsembleSV/scripts
cp soft/EnsembleSV/data.yaml .
cp soft/EnsembleSV/sv_tools.yaml .
````

Now update the copied `data.yaml` and `sv_tools.yaml` files with the experiment-specific information. 
On detailed instruction for updating `data.yaml` file, please, refer to respective data [docs](./docs/data.md).
On detailed instruction for updating `sv_tools.yaml` file, please, refer to respective tools [docs](./docs/sv_tools.md).

Running *EnsembleSV* can be accomplished via *snakemake* simple command:
````
snakemake -s merge_svs.snakefile
````

If you want to separate method-specific SV calling and subsequent merging, you can do so as follows:
````
snakemake -s call_svs.snakefile
snakemake -s merge_svs.snakefile
```` 

Running method-specific SV calling can be achieved via:
````
snakemake -s call_svs_*method*.snakefile
````

Note that for every data type (short Illumina, linked, and long reads) only SV inference methods specified in the `tools_enabled_methods` section in the `sv_tools.yaml` file.  

Contribution
=
If you wish to contribute to the EnsembleSV project, please, contact Sergey Aganezov via email [sergeyaganezovjr(at)gmail.com](mailto:segreyaganezovjr@gmail.com) or submit a pull request with suggested additions/changes.

Issues
=
If you identify any issues and/or bugs with the *EnsembleSV* pipeline, or want to suggest an enhancement to it, please, use the repository-associated [issue tracker](https://github.com/aganezov/EnsembleSV/issues).  

Citation
=
If you use *EnsembleSV* in your research, please cite the following manuscript [TBA]().