## Additional flags to pass to aclocal when it is invoked automatically at
## make time. The ${ACLOCAL_FLAGS} variable is picked up from the environment
## to provide a way for the user to supply additional arguments.
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

## Paths to include
TABIX_VER = 0.2.5
TABIX = tabix-$(TABIX_VER)
AM_CPPFLAGS = -Wall -Wno-unused-local-typedefs -Wno-enum-compare -fpic -O2 #-I../include/ #-I$(TABIX)/

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

bin_PROGRAMS = vcfast bgzip tabix tinyfast

noinst_LIBRARIES = libtabix.a

## vcfast files
vcfast_SOURCES = vcfast.cpp Error.cpp PhredHelper.cpp params.cpp codonHelper.cpp genomePosition.cpp IO.cpp
vcfast_LDADD = libtabix.a

tinyfast_SOURCES = tinyfast.cpp Error.cpp params.cpp
tinyfast_LDADD = libtabix.a

bgzip_SOURCES = bgzip.c bedidx.c bgzf.c knetfile.c kstring.c

tabix_SOURCES = tabix.c bedidx.c bgzf.c knetfile.c kstring.c index.c

libtabix_a_SOURCES = bedidx.c bgzf.c knetfile.c kstring.c index.c
