## Additional flags to pass to aclocal when it is invoked automatically at
## make time. The ${ACLOCAL_FLAGS} variable is picked up from the environment
## to provide a way for the user to supply additional arguments.
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

#SUBDIRS = htslib

## Paths to include
#TABIX_VER = 0.2.5
#TABIX = tabix-$(TABIX_VER)
AM_CPPFLAGS = -fopenmp -I ../htslib/ -I ../htslib/htslib -pipe -D__STDC_LIMIT_MACROS -Wall -Wno-unused-local-typedefs -Wno-enum-compare -fpic -O2 #-I../include/ #-I$(TABIX)/
#AM_CPPFLAGS = -I ../htslib/ -I ../htslib/htslib -pipe -D__STDC_LIMIT_MACROS -Wall -Wno-unused-local-typedefs -Wno-enum-compare -fpic -g #-O2 #-I../include/ #-I$(TABIX)/

noinst_HEADERS = hDist.h \
	bam_endian.h \
	boolParser.h \
	cdsStat.h ncStat.h \
	TypeConversion.h \
	Constant.h GenomeSequence.h \
	kseq.h pException.h tabix.h \
	Error.h \
	ksort.h pFile.h \
	genomeLoci.h genomePosition.h IO.h kstring.h \
	PhredHelper.h bgzf.h \
	codonHelper.h fVcf.h \
	genomeScore.h khash.h params.h \
	boolExprParser.h \
	knetfile.h \
	wFile.h 

#bin_PROGRAMS = vcfast bgzip tabix tinyfast
bin_PROGRAMS = vcfast tinyfast bcfast cramore

noinst_LIBRARIES = libtabix.a

## vcfast files
cramore_SOURCES = cramore.cpp \
	Error.cpp PhredHelper.cpp params.cpp \
	dropseq.cpp \
	hts_utils.cpp bcf_ordered_reader.cpp interval.cpp interval_tree.cpp utils.cpp genome_interval.cpp log_tool.cpp estimator.cpp bcf_ordered_writer.cpp \
	bam_ordered_reader.cpp \
	nuclear_pedigree.cpp \
	genotype_concordance.cpp \
	filter.cpp \
	cmd_vcf_mendel_dup_conc.cpp cmd_vcf_sample_summary.cpp cmd_vcf_squeeze.cpp \
	cmd_cram_simul_contam.cpp cmd_cram_verify_pair_id.cpp cmd_cram_sparse_genotype.cpp cmd_cram_flagstat.cpp \
	cmd_cram_procap_detect.cpp \
	cmd_sc_multinom_em.cpp cmd_sc_map_stamps.cpp cmd_sc_kallisto_count.cpp 
cramore_LDADD = ../htslib/libhts.a -lpthread
#cramore_LDFLAGS = -fopenmp

## vcfast files
bcfast_SOURCES = bcfast.cpp Error.cpp PhredHelper.cpp params.cpp bcf_ordered_reader.cpp interval.cpp interval_tree.cpp utils.cpp hts_utils.cpp genome_interval.cpp
bcfast_LDADD = ../htslib/libhts.a -lpthread

## vcfast files
vcfast_SOURCES = vcfast.cpp Error.cpp PhredHelper.cpp params.cpp codonHelper.cpp genomePosition.cpp IO.cpp
vcfast_LDADD = libtabix.a

tinyfast_SOURCES = tinyfast.cpp Error.cpp params.cpp
tinyfast_LDADD = libtabix.a

#bgzip_SOURCES = bgzip.c bedidx.c bgzf.c knetfile.c kstring.c

#tabix_SOURCES = tabix.c bedidx.c bgzf.c knetfile.c kstring.c index.c

libtabix_a_SOURCES = bedidx.c bgzf.c knetfile.c kstring.c index.c
