#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([apigenome], [0.1.0], [hmkang@umich.edu])

AC_CONFIG_AUX_DIR(aux-dist)

#AM_INIT_AUTOMAKE
AM_INIT_AUTOMAKE([subdir-objects])
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_MACRO_DIR([m4])

#AM_PROG_LIBTOOL
#AC_PROG_LIBTOOL
#AC_ENABLE_SHARED
#AC_DISABLE_STATIC

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_F77

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([limits.h stddef.h stdint.h stdlib.h string.h unistd.h])

#AC_SUBST(CPPFLAGS, "$CPPFLAGS -I ${abs_top_builddir}/../htslib")
#AC_SUBST(LDFLAGS, "$LDFLAGS -L ${abs_top_builddir}/../htslib -lhts")

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_C_RESTRICT
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_HEADER_MAJOR
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor gethostname memset pow select sqrt strchr strstr strtol])
AC_CHECK_HEADERS([zlib.h], [], [AC_MSG_ERROR([zlib.h is not found in the inclusion path. In Ubuntu, it can be installed by running 'sudo apt-get install libz-dev'])])
AC_CHECK_LIB([z], [gzopen], [], [AC_MSG_ERROR([libz.{so,a} is not found in the library path. In Ubuntu, it can be installed by running 'sudo apt-get install libz-dev'])])
AC_CHECK_LIB([bz2], [BZ2_bzopen], [], [AC_MSG_ERROR([libbz2.{so,a} was not found in the library path. In Ubuntu, it can be installed by running 'sudo apt-get install bzip2-devel'])])
AC_CHECK_LIB([crypto], [EVP_EncryptInit], [], [AC_MSG_ERROR([openssl crypto library is not found in the library path.  In Ubuntu, it can be installed by running 'sudo apt-get install openssl libssl-dev'])])
AC_CHECK_LIB([ssl], [SSL_CTX_new], [], [AC_MSG_ERROR([SSL library is not found in the library path. In Ubuntu, it can be installed by running 'sudo apt-get install openssl libssl-dev'])])
AC_CHECK_LIB([curl],[curl_easy_init],[], [AC_MSG_ERROR([libcurl is not found in the library path. In Ubuntu, it can be installed by running 'sudo pt-get install libcurl4-openssl-dev'])])
AC_CHECK_LIB([lzma], [lzma_code], [], [AC_MSG_ERROR([liblzma is not found in the library path. In Ubuntu, it can be installed by running 'sudo apt-get install lzma'])])

#AC_CHECK_HEADERS([omp.h], [], [AC_MSG_ERROR([omp.h is not found in the inclusion path. In Ubuntu, it can be installed by running 'sudo apt-get install libomp-dev'])])
#AC_CHECK_LIB([gomp],[omp_get_num_threads], [], [AC_MSG_ERROR([libgomp is not found in the library path. In Ubuntu, it can be installed by running 'sudo apt-get install libomp-dev'])])

AC_CHECK_HEADERS([math.h], [], [AC_MSG_ERROR([math.h is not found in the inclusion path. In Ubuntu, it can be installed by running 'sudo apt-get install libc6-dev'])])
AC_CHECK_LIB([m],[erf],[AC_DEFINE([HAVE_ERF],[1],[libm includes erf])])
AC_CHECK_LIB([m],[erfc],[AC_DEFINE([HAVE_ERFC],[1],[libm includes erfc])])

#AC_CHECK_HEADERS([htslib/bgzf.h], [], [AC_MSG_ERROR([htslib/hts.h is not found in the inclusion path (typically under the parent directory. Visit http://github.com/samtools/htslib for installation and setup])])
#AC_CHECK_LIB([hts], [hts_open], [], [AC_MSG_ERROR([htslib is not found in the library path (typically under the parent directory). Visit http://github.com/samtools/htslib for installation and setup])])
#AC_CHECK_LIB(R, R_alloc, [], [AC_MSG_ERROR([libR.{so,a} was not found. Please install it at http://www.r-project.org/ first])])

LT_PREREQ([2.2])
LT_INIT([dlopen shared])

#perl 
AC_CHECK_PROG(PERL_CHECK,perl,yes)
if test x"$PERL_CHECK" != x"yes" ; then
   AC_MSG_ERROR([Please install perl before installing.])
fi

#gs
#AC_CHECK_PROG(GS_CHECK,gs,yes)
#if test x"$GS_CHECK" != x"yes" ; then
#   AC_MSG_ERROR([Please install gs before installing.])
#fi

#Rscript
AC_CHECK_PROG(RSCRIPT_CHECK,Rscript,yes)
if test x"$RSCRIPT_CHECK" != x"yes" ; then
   AC_MSG_ERROR([Please install latest version of R before installing.])
fi

#R
AC_CHECK_PROG(R_CHECK,R,yes)
if test x"$R_CHECK" != x"yes" ; then
   AC_MSG_ERROR([Please install latest version of R before installing.])
fi

#grep
AC_CHECK_PROG(GREP_CHECK,grep,yes)
if test x"$GREP_CHECK" != x"yes" ; then
   AC_MSG_ERROR([Please install grep before installing.])
fi

#mv
AC_CHECK_PROG(MV_CHECK,mv,yes)
if test x"$MV_CHECK" != x"yes" ; then
   AC_MSG_ERROR([Please install mv before installing.])
fi

#rm
AC_CHECK_PROG(RM_CHECK,rm,yes)
if test x"$RM_CHECK" != x"yes" ; then
   AC_MSG_ERROR([Please install rm before installing.])
fi

#make
AC_CHECK_PROG(MAKE_CHECK,make,yes)
if test x"$MAKE_CHECK" != x"yes" ; then
   AC_MSG_ERROR([Please install make before installing.])
fi

#cat
AC_CHECK_PROG(CAT_CHECK,cat,yes)
if test x"$CAT_CHECK" != x"yes" ; then
   AC_MSG_ERROR([Please install cat from GNU coreutils at http://ftp.gnu.org/gnu/coreutils/ before installing.])
fi

#cut
AC_CHECK_PROG(CUT_CHECK,cut,yes)
if test x"$CUT_CHECK" != x"yes" ; then
   AC_MSG_ERROR([Please install cut from GNU coreutils at http://ftp.gnu.org/gnu/coreutils/ before installing.])
fi

#dirname
AC_CHECK_PROG(DIRNAME_CHECK,dirname,yes)
if test x"$DIRNAME_CHECK" != x"yes" ; then
   AC_MSG_ERROR([Please install dirname from GNU coreutils at http://ftp.gnu.org/gnu/coreutils/ before installing.])
fi

#head
AC_CHECK_PROG(HEAD_CHECK,head,yes)
if test x"$HEAD_CHECK" != x"yes" ; then
   AC_MSG_ERROR([Please install head from GNU coreutils at http://ftp.gnu.org/gnu/coreutils/ before installing.])
fi

#mkdir
AC_CHECK_PROG(MKDIR_CHECK,mkdir,yes)
if test x"$MKDIR_CHECK" != x"yes" ; then
   AC_MSG_ERROR([Please install mkdir from GNU coreutils at http://ftp.gnu.org/gnu/coreutils/ before installing.])
fi

#sort
AC_CHECK_PROG(SORT_CHECK,sort,yes)
if test x"$SORT_CHECK" != x"yes" ; then
   AC_MSG_ERROR([Please install sort from GNU coreutils at http://ftp.gnu.org/gnu/coreutils/ before installing.])
fi

#zcat
AC_CHECK_PROG(ZCAT_CHECK,zcat,yes)
if test x"$ZCAT_CHECK" != x"yes" ; then
   AC_MSG_ERROR([Please install zcat from zlib package at http://www.zlib.net/ before installing.])
fi

#groff
#AC_CHECK_PROG(GROFF_CHECK,groff,yes)
#if test x"$GROFF_CHECK" != x"yes" ; then
#   AC_MSG_ERROR([Please install groff at ftp://ftp.gnu.org/gnu/groff/ before installing.])
#fi

#pfbtops
#AC_CHECK_PROG(PFBTOPS_CHECK,pfbtops,yes)
#if test x"$PFBTOPS_CHECK" != x"yes" ; then
#   AC_MSG_ERROR([Please install pfbtops from groff package at ftp://ftp.gnu.org/gnu/groff/ before installing.])
#fi

#gnuplot
#AC_CHECK_PROG(GNUPLOT_CHECK,gnuplot,yes)
#if test x"$GNUPLOT_CHECK" != x"yes" ; then
#   AC_MSG_ERROR([Please install gnuplot at http://sourceforge.net/projects/gnuplot/files/ before installing.])
#fi

AC_CONFIG_FILES([Makefile libnmf/Makefile src/Makefile scripts/Makefile data/Makefile])

#AC_OUTPUT(Makefile libnmf/Makefile src/Makefile scripts/Makefile data/Makefile)
AC_OUTPUT
