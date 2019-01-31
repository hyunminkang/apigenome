# APIGenome - Big data genomics analysis libraries & tools

### Introduction

**_APIGenome_** consists of a collection of perl libraries and command-line utilities for big data genomics analysis. This repository started as a personal repository of small genomic analysis tools developed by [Hyun Min Kang](https://sph.umich.edu/faculty-profiles/kang-hyunmin.html), but some of the utilities developed in this repository may be exposed more widely.

**_APIGenome_** provides a useful perl software development, in regards to handling command-line arguments, automated self-documentation, and useful APIs for genomic analysis.

Because **_APIGenome_** contains many _"in-progress"_ software tools as a _sandbox_ (as the name of default repository indicates), it contains many software tools that has not been fully described or documented. As APIGenome will be continuously updated with new tools, it will be always **_UNDER CONSTRUCTION_**. Many preliminary utilities will be partially documented, and may contain bugs. So use the software tools at your own risk.

Note that some of the software tools in this repository **_MAY MIGRATE OUT TO A STANDALONE REPOSITORY_** if the tool receives a wide attention enough to have their own brand name. 


### Installing APIGenome

#### Requirement before installation

Currently, APIGenome installation was tested in Ubuntu and Mac OS X. If you find installation problems in other OS, please let the authors know.

First, you need standard UNIX tools including `grep, mv, rm, make, cat, cut, dirname, head, mkdir, sort, zcat` installed.

Next, you will need to have a number of tools installed, including `autoconf, automake, libtool, perl, R, Rscript`.


#### Cloning from github repository

At the parent directory of `htslib/`, you can clone the current snapshot of this repository to install as well. 

<pre>
$ git clone https://github.com/hyunminkang/apigenome.git
$ cd apigenome/
$ autoreconf -vfi
$ ./configure --prefix [/path/to/install]
$ make
$ make install </pre>

It is recommended to specify `--prefix` explictly as the installation without `--prefix` argument has not been extensively tested.

### How to use APIGenome utilities

APIGenome contains a list of many self-documented command line utilities. To understand how to use each of them, you can run each utility with -man or -help option to see the command line usages.

<pre>
$ [path/to/apigenome]/bin/[utility-name] -man 
$ [path/to/apigenome]/bin/[utility-name] -help
</pre>

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
