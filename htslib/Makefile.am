## Additional flags to pass to aclocal when it is invoked automatically at
## make time. The ${ACLOCAL_FLAGS} variable is picked up from the environment
## to provide a way for the user to supply additional arguments.
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

## Paths to include
AM_CPPFLAGS = -pthread -Wall -Wno-unused-local-typedefs -Wno-enum-compare -fpic -O2

noinst_HEADERS = \
	cram/cram.h \
	cram/cram_codecs.h \
	cram/cram_decode.h \
	cram/cram_encode.h \
	cram/cram_index.h \
	cram/cram_io.h \
	cram/cram_samtools.h \
	cram/cram_stats.h \
	cram/cram_structs.h \
	cram/mFILE.h \
	cram/misc.h \
	cram/open_trace_file.h \
	cram/os.h \
	cram/pooled_alloc.h \
	cram/rANS_byte.h \
	cram/rANS_static.h \
	cram/sam_header.h \
	cram/string_alloc.h \
	cram/vlen.h \
	cram/zfio.h \
	hfile_internal.h \
	hts_internal.h \
	htslib/bgzf.h \
	htslib/cram.h \
	htslib/faidx.h \
	htslib/hfile.h \
	htslib/hts.h \
	htslib/hts_defs.h \
	htslib/hts_endian.h \
	htslib/kbitset.h \
	htslib/kfunc.h \
	htslib/khash.h \
	htslib/khash_str2int.h \
	htslib/klist.h \
	htslib/knetfile.h \
	htslib/kseq.h \
	htslib/ksort.h \
	htslib/kstring.h \
	htslib/regidx.h \
	htslib/sam.h \
	htslib/synced_bcf_reader.h \
	thread_pool_internal.h \
	htslib/thread_pool.h \
	htslib/tbx.h \
	htslib/vcf.h \
	htslib/vcf_sweep.h \
	htslib/vcfutils.h \
	config.h \
	version.h 

#bin_PROGRAMS = vcfast bgzip tabix tinyfast
bin_PROGRAMS = bgzip tabix

noinst_LIBRARIES = libhts.a

bgzip_SOURCES = bgzip.c
bgzip_LDADD = libhts.a -lpthread

tabix_SOURCES = tabix.c
tabix_LDADD = libhts.a -lpthread

libhts_a_SOURCES = \
	bgzf.c \
	errmod.c \
	faidx.c \
	hfile.c \
	hfile_gcs.c \
	hfile_libcurl.c \
	hfile_net.c \
	hfile_s3.c \
	hts.c \
	kfunc.c \
	knetfile.c \
	kstring.c \
	md5.c \
	multipart.c \
	plugin.c \
	probaln.c \
	realn.c \
	regidx.c \
	sam.c \
	synced_bcf_reader.c \
	vcf_sweep.c \
	tbx.c \
	textutils.c \
	thread_pool.c \
	vcf.c \
	vcfutils.c \
	cram/cram_codecs.c \
	cram/cram_decode.c \
	cram/cram_encode.c \
	cram/cram_external.c \
	cram/cram_index.c \
	cram/cram_io.c \
	cram/cram_samtools.c \
	cram/cram_stats.c \
	cram/files.c \
	cram/mFILE.c \
	cram/open_trace_file.c \
	cram/pooled_alloc.c \
	cram/rANS_static.c \
	cram/sam_header.c \
	cram/string_alloc.c \
	cram/vlen.c \
	cram/zfio.c
