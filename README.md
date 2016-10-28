# APIGenome - A public library of big data genomics analysis tools

### What is APIGenome?

APIGenome is a collection of libraries and command-line utilities for big data genomics analysis. It contains a mixture of mature and nascent software tools actively used by research groups. APIGenome provides a useful software framework for developers to self-document their tools. Currently, it is developed and maintained by Hyun Min Kang's group.

As APIGenome will be continuously updated with new tools, it will be always **UNDER CONSTRUCTION**. Many preliminary utilities will be partially documented, and may contain bugs.

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

The list of currently provided utilities include:

| Category | Utility Name  | Description |
| :-----: |:--------------:| :----------------------- |
| 
| Sequence Reads | demux-fastq | Demultiplex barcoded FASTQ for single-ended sequence data in DropSeq format |
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
|               | vcf-liftover | Software for lifting over VCF files |
|               | vcf-delta-svm | DeltaSVM Method (Lee et al. Nat Genet 2015) implemented for VCF files |
|               | vcf-summary | Produce basic summary of VCF files such as Ts/Tv, %dbSNP |
| Expression | cel-extract-intensity | Extract probe-level intensity information from CEL files |
| PROcap/PROseq | pileup-pro | Produce pileups for PROseq/PROcap data |
|               | bed-tss-match | Identify transcription start sites from pileups |
|                | rev-trim | Reverse complement and trim sequence reads |
| DropSeq | align-dropseq | All-in-one alignment of DropSeq sequence reads |
| Parallelization | run-make | Fault-tolerant parallization utility based on Makefile |
