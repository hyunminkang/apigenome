# APIGenome - A public library of big data genomics analysis tools

### What is APIGenome?

APIGenome is a collection of libraries and command-line utilities for big data genomics analysis. It contains a mixture of mature and nascent software tools actively used by research groups. APIGenome provides a useful software framework for developers to self-document their tools. Currently, it is developed and maintained by Hyun Min Kang's group.

As APIGenome will be continuously updated with new tools, it will be always UNDER CONSTRUCTION. Many nascent utilities will be partially documented, and may contain bugs.

### Installing APIGenome

#### Downloading the from latest release

To install APIGenome, you can download the [Latest Release of APIGenome](https://github.com/hyunminkang/apigenome/raw/master/releases/apigenome-latest.tar.gz) by clicking the link or using command line.
<pre>
$ wget https://github.com/hyunminkang/apigenome/raw/master/releases/apigenome-latest.tar.gz </pre>

Note that the latest release is a stable snapshot of a previous version of APIGenome, but may not contain the more recent update since the last freeze.

To install APIGenome, uncompress the tarball, configure the install directory, build, and install.

<pre>
$ tar xzvf apigenome-latest.tar.gz
$ cd apigenome-[version]/
$ ./configure --prefix [/path/to/install]
$ make
$ make install </pre>

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

### List of stable APIGenome utilities

APIGenome contains many utilities that are not finished and maybe for internal use. You may use them at your own risk too, but the software tools listed below should be mature enough to get assistance from the developers when help is needed.

List of currently provided utilities:

| Category | Utility Name  | Description |
| :-----: |:--------------:| :----------------------- |
| 
| Sequence Reads | demux-fastq | Demultiplex barcoded FASTQ for single-ended sequence data |
|                | bam-quick-peek-batch | Produce simple summary statistics for a list of BAM files |
|                | now-seq-batch | Produce a summary of QC metrics from outputs of GotCloud alignment pipeline |
| Variant Calls | vcfast | Fast command line utility for processing VCF files (in C++) |
|               | vcf-summary | Produce basic summary of VCF files such as Ts/Tv, %dbSNP |
|               | vcf-milk-filter | Mendelian-inheritance and likelihood based variant filtering software |
|               | vcf-add-rsid | Add rsIDs to VCF files |
|               | vcf-lookup-rsid | Lookup a variant from VCF based on rsIDs |
|               | vcf-resolve-chrX-hets | Resolve heterozygous genotyes from chrX based on likelihoods |
|               | vcf-issac-sumary | Sumamrize VCF files produced by Illumina's iSSAC pipeline |
|               | draw-afs | Draw allele frequency spectrum (AFS) from VCF sites |
| Expression | cel-extract-intensity | Extract probe-level intensity information from CEL files |
| PROcap/PROseq | pileup-pro | Produce pileups for PROseq/PROcap data |
|               | bed-tss-match | Identify transcription start sites from pileups |
|                | rev-trim | Reverse complement and trim sequence reads |
| DropSeq | align-dropseq | All-in-one alignment of DropSeq sequence reads |
| Parallelization | run-make | Fault-tolerant parallization utility based on Makefile |
