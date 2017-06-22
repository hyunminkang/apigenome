# APIGenome - Big data genomics analysis libraries & tools

### NOTE TO USERS OF demuxlet
If you are interested in installing **_demuxlet_** software tool, please visit the standalone release http://github.com/statgen/demuxlet . This repository can be under construction often, so may not build smoothly. 

### Introduction

**_APIGenome_** consists of a collection of perl/C++ libraries and command-line utilities for big data genomics analysis. This repository started as a personal repository of small genomic analysis tools developed by [Hyun Min Kang](https://sph.umich.edu/faculty-profiles/kang-hyunmin.html), but some of the utilities developed in this repository may be exposed more widely.

**_APIGenome_** provides a useful perl/C++ software development, in regards to handling command-line arguments, automated self-documentation, and useful APIs for genomic analysis.

Because **_APIGenome_** contains many _"in-progress"_ software tools as a _sandbox_ (as the name of default repository indicates), it contains many software tools that has not been fully described or documented. As APIGenome will be continuously updated with new tools, it will be always **_UNDER CONSTRUCTION_**. Many preliminary utilities will be partially documented, and may contain bugs. So use the software tools at your own risk.

Note that some of the software tools in this repository **_MAY MIGRATE OUT TO A STANDALONE REPOSITORY_** if the tool receives a wide attention enough to have their own brand name. 


### Installing APIGenome

#### Cloning from github repository

You can clone the current snapshot of this repository to install as well

<pre>
$ git clone https://github.com/hyunminkang/apigenome.git
$ cd apigenome/
$ ./configure --prefix [/path/to/install]
$ make
$ make install </pre>

### How to use APIGenome utilities

APIGenome contains a list of many self-documented command line utilities. To understand how to use each of them, you can run each utility with -man or -help option to see the command line usages.

<pre>
$ [path/to/apigenome]/bin/[utility-name] -man </pre>

### List of available APIGenome utilities

The complete list of APIGenome Utility can be found at 
<pre>$(INSTALL_DIR)/bin/</pre> 
Some of these utilities that are not finished and maybe under development for internal use. You may use them at your own risk, but the software tools listed below should be relatively more mature enough to get assistance from the developers when help is needed.

Each software tool is self-documented, and you can see the detailed document by typing
<pre>$(INSTALL_DIR)/bin/$(PROGRAM_NAME) -man</pre>

For C++ program, such as <code>cramore</code> contains multiple programs inside it, and the documentation can be found by typing
<pre>$(INSTALL_DIR)/bin/cramore $(COMMAND_NAME) -help</pre>

The development status of the software are classified into four stages:
* alpha : The software tool is in early stage of development, or it is not yet ready for sharing. It may be buggy and the documentation is not comprehensive. Developer may not promptly respond to bug report or questions. Use at your own (high) risk. 
* beta : The software tool is in a relatively good shape to share with others. It is either (a) in the stage of testing or (b) considered as a narrowly-shared (in-house) tools. It is self-documented, but it is unlikely to be comprehensively documented and interfaced because the software does not have wide attention, or it is in the test stage of development. For the software tool that is expected to move to release status, the documentation and interface will continute to be improved. Bugs may exist, but bug reports will receive attention and appreciation by authors. 
* release : The software tool is released, and possibly already published. Polishing documentation, bug reports will be much appreciated. Questions regarding the released software tools will be more promptly answered.
* migrating : The software tool received a wide attention enough to be justified to be migrated and branded on its own. During the migration period, there can be two separate repository of the tool, but the new repository that will be migrated into may be more up-to-date, and the version in this repository is likely to stay in a more stable version.

Below is the list of **_APIGenome_** utilities, except for `cramore` software tool. Note that the list is only occasionally updated, so may not reflect recent changes. See above if you need a more detailed documentation

| Category | Utility Name  | Development Status | Brief Description |
| :------ |:--------------:|:--------------:|:----------------------- |
| 
| Sequence Reads | `align-dropseq` | alpha | All-in-one alignment of DropSeq sequence reads |
|                | `align-pro` | alpha | Align PROcap sequence data |
|                | `bam-quick-peek-batch` | alpha | Produce simple summary statistics for a list of BAM files |
|                | `dropseq-resolve-barcode-trimming` | alpha | Resolve barcode trimming issue for Dropseq |
|                | `dropseq-kallisto` | alpha | Sequence alignment of Dropseq data using kallisto software tool |
|                | `demux-fastq` | alpha | Demultiplex barcoded FASTQ for single-ended sequence data in DropSeq format |
|                | `now-seq-batch` | alpha | Produce a summary of QC metrics from outputs of GotCloud alignment pipeline |
|                | `pileup-pro` | alpha | Produce pileups for PROseq/PROcap data |
|                | `rev-trim` | alpha | Reverse complement and trim sequence reads |
| Variant Calls | `bed-diff` | beta | Compare genotype concordance and discordance between callset |
|               | `draw-afs` | alpha | Draw allele frequency spectrum (AFS) from VCF sites |
|               | `inspect-sv` | alpha | Examine the reads within a structural variant |
|               | `vcf-lookup-rsid` | beta | Lookup a variant from VCF based on rsIDs |
|               | `vcf-issac-sumary` | alpha | Sumamrize VCF files produced by Illumina's iSSAC pipeline |
|               | `vcf-add-rsid` | beta | Add rsIDs to VCF files |
|               | `vcf-delta-svm` | beta | DeltaSVM Method (Lee et al. Nat Genet 2015) implemented for VCF files |
|               | `vcf-extract-gt-only` | alpha | Extract GT field only from a VCF |
|               | `vcf-f2-sharing` | alpha | Extract GT field only from a VCF |
|               | `vcf-liftover` | beta | Software for lifting over VCF files |
|               | `vcf-lookup-rsid` | beta | Lookup rsID from a VCF |
|               | `vcf-milk-filter` | alpha | Mendelian-inheritance and likelihood based variant filtering software |
|               | `vcf-resolve-chrX-hets` | alpha | Resolve heterozygous genotyes from chrX based on likelihoods |
|               | `vcf-summary` | beta | Produce basic summary of VCF files such as Ts/Tv, %dbSNP |
|               | `vcf-summary-merge` | alpha | Merge multiple `vcf-summary` output files |
|               | `vcf-summary-v2` | alpha | A newer version of `vcf-summary` that includes indels |
|               | `vcfast` | beta | Fast command line utility for processing VCF files (in C++) |
| Expression Array | `cel-extract-intensity` | alpha | Extract probe-level intensity information from CEL files |
| Other Genomics Tool | `bed-tss-match` | alpha | Identify transcription start sites from PROcap pileups |
|                | `count-dropseq` | alpha | Produce digital expression matrix from Dropseq |
|                 | `gene-conv-name` | alpha | convert between gene names and symbols |
|                 | `merge-tsv-generic` | beta | Sum the values across multiple tab-delimited files formatted in a same way |
|                 | `tsv-join` | beta | Join multiple tab-delimited files pivoted by a shared column |
| Parallelization | `run-make` | beta | Fault-tolerant parallization utility based on Makefile |


Below is the list of C++ utilities supported by `cramore`. Note that you need to run `$(INSTALL_DIR)/bin/cramore UTILITY_NAME` to run them. 

| Category | Utility Name  | Development Status | Brief Description |
| :------ |:--------------:|:--------------:|:----------------------- |
| 
| BED           | `bed-delta-svm-train` | alpha | Train delta-SVM using lsgkm |
|               | `bed-matched-shuffle` | alpha | Shuffle genomic intervals |
| SAM/BAM/CRAM | `cram-simul-contam` | alpha | Simulate contaminated reads |
|                | `cram-context-indel-analysis` | alpha | Context-specific indel analysis |
|                | `cram-flagstat-all` | beta | Produce a comprehensive flagstat summary |
|                | `cram-procap-detect` | alpha | TSS detection of PROcap data |
|                | `cram-sparse-genotype` | alpha | Sparse genotyping from SAM and VCF |
|                | `demuxlet`          | beta | Sample demultiplexing of dsc-RNAseq reads |
|                | `verify-pair-id`          | alpha | Initial version of `demuxlet` |
| BCF/VCF        | `vcf-delta-svm`  | beta | Newer version of delta-svm compatible to lsgkm |
|                | `vcf-extract`  | alpha | Extract specific variants from VCF files |
|                | `vcf-infer-ancestry`  | alpha | Infer ancestry from VCF |
|                | `vcf-infer-isaf`  | alpha | Infer individual-specific allele frequencies from VCF |
|                | `vcf-mendel-dup-conc`  | beta | Calculate Mendelian/duplicate genotype concordance |
|                | `vcf-sample-summary`  | beta | Produce per-sample summary from BCF/VCF |
|                | `vcf-squeeze`  | beta | Extract partial information from BCF/VCF |
|                | `vcf-svd`      | alpha | Perform singular value decomposition from VCF |
|                | `vcf-update-sites`      | alpha | Update VCF site information using another VCF file |
| TXT            | `multinom-em` | alpha | Perform multinomial EM on digital expression matrix
