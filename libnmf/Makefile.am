## Additional flags to pass to aclocal when it is invoked automatically at
## make time. The ${ACLOCAL_FLAGS} variable is picked up from the environment
## to provide a way for the user to supply additional arguments.
ACLOCAL_AMFLAGS = -Iinclude/ m4 ${ACLOCAL_FLAGS}

## Paths to include
AM_CPPFLAGS = -Iinclude/ -Wall -Wextra -Wno-unused-local-typedefs -Wno-enum-compare -fpic -static -O2 

noinst_HEADERS = \
	include/blaslapack.h \
	include/calculatemaxchange.h \
	include/calculatemaxchange_singleprec.h \
	include/calculatenorm.h \
	include/calculatenorm_singleprec.h \
	include/calculatesvd.h \
	include/calculatesvd_singleprec.h \
	include/checkarguments.h \
	include/checkarguments_singleprec.h \
	include/checkmatrices.h \
	include/checkmatrices_singleprec.h \
	include/common.h \
	include/cpmidx_double.h \
	include/cpmidx_float.h \
	include/generatematrix.h \
	include/generatematrix_singleprec.h \
	include/getattrvalvsclassoccurence.h \
	include/getattrvalvsclassoccurence_singleprec.h \
	include/getentropy.h \
	include/getentropy_singleprec.h \
	include/getinfogainforattribute.h \
	include/getinfogainforattribute_singleprec.h \
	include/infogain.h \
	include/infogain_openmp.h \
	include/infogain_openmp_singleprec.h \
	include/infogain_singleprec.h \
	include/loadmatrix.h \
	include/loadmatrix_singleprec.h \
	include/nmf_als.h \
	include/nmf_alspg.h \
	include/nmf_alspg_singleprec.h \
	include/nmf_als_singleprec.h \
	include/nmf_bayes.h \
	include/nmf_bayes_singleprec.h \
	include/nmfdriver.h \
	include/nmfdriver_singleprec.h \
	include/nmf_fast.h \
	include/nmf_fast_singleprec.h \
	include/nmf_mu.h \
	include/nmf_mu_singleprec.h \
	include/nmf_neals.h \
	include/nmf_neals_singleprec.h \
	include/nmf_pg.h \
	include/nmf_pg_singleprec.h \
	include/outputtiming.h \
	include/pg_subprob_h.h \
	include/pg_subprob_h_singleprec.h \
	include/pg_subprob_w.h \
	include/pg_subprob_w_singleprec.h \
	include/randnumber.h \
	include/randnumber_singleprec.h \
	include/setdefaultopts.h \
	include/setdefaultopts_singleprec.h \
	include/storematrix.h \
	include/storematrix_singleprec.h

bin_PROGRAMS = nmf

noinst_LIBRARIES = libnmf.a libarpack.a

nmf_SOURCES = nmf_main.c
nmf_LDADD = libnmf.a libarpack.a -llapack -lblas -lm
nmf_LDFLAGS = -fopenmp

libarpack_a_SOURCES = \
	idamax.f \
	dasum.f \
	daxpy.f \
	dcopy.f \
	ddot.f \
	dnrm2.f \
	drot.f \
	drotg.f \
	dscal.f \
	dswap.f \
	dgemv.f \
	dger.f \
	dsymv.f \
	dsyr2.f \
	dtrmm.f \
	dgbmv.f \
	dgeqr2.f \
	dlabad.f \
	dlacon.f \
	dlacpy.f \
	dladiv.f \
	dlae2.f \
	dlaev2.f \
	dlaexc.f \
	dlagtm.f \
	dlahqr.f \
	dlaln2.f \
	dlamch.f \
	dlange.f \
	dlanhs.f \
	dlanst.f \
	dlanv2.f \
	dlaptm.f \
	dlapy2.f \
	dlapy3.f \
	dlaran.f \
	dlarf.f \
	dlarfg.f \
	dlarfx.f \
	dlarnd.f \
	dlarnv.f \
	dlartg.f \
	dlaruv.f \
	dlascl.f \
	dlaset.f \
	dlasr.f \
	dlasrt.f \
	dlassq.f \
	dlasy2.f \
	dorm2r.f \
	dsteqr.f \
	dtrevc.f \
	dtrexc.f \
	dtrsen.f \
	dtrsyl.f \
	ilaenv.f \
	lsame.f \
	lsamen.f \
	xerbla.f \
	xlaenv.f \
	dvout.f \
	dmout.f \
	icnteq.f \
	icopy.f \
	iset.f \
	iswap.f \
	ivout.f \
	second.f \
	dgetv0.f \
	dlaqrb.f \
	dstqrb.f \
	dsortc.f \
	dsortr.f \
	dstatn.f \
	dstats.f \
	dnaitr.f \
	dnapps.f \
	dnaup2.f \
	dnaupd.f \
	dnconv.f \
	dneigh.f \
	dngets.f \
	dsaitr.f \
	dsapps.f \
	dsaup2.f \
	dsaupd.f \
	dsconv.f \
	dseigt.f \
	dsgets.f \
	dneupd.f \
	dseupd.f \
	dsesrt.f

libnmf_a_SOURCES = \
	checkarguments.c \
	checkarguments_singleprec.c \
	infogain_openmp.c \
	infogain_openmp_singleprec.c \
	infogain.c \
	infogain_singleprec.c \
	getentropy.c \
	getentropy_singleprec.c \
	getinfogainforattribute.c \
	getinfogainforattribute_singleprec.c \
	getattrvalvsclassoccurence.c \
	getattrvalvsclassoccurence_singleprec.c \
	cpmidx_double.c \
	cpmidx_float.c \
	loadmatrix.c \
	loadmatrix_singleprec.c \
	nmfdriver.c \
	nmfdriver_singleprec.c \
	nmf_pg.c \
	nmf_pg_singleprec.c \
	nmf_fast.c \
	nmf_fast_singleprec.c \
	nmf_bayes.c \
	nmf_bayes_singleprec.c \
	randnumber.c \
	randnumber_singleprec.c \
        calculatemaxchange.c \
	calculatemaxchange_singleprec.c \
	checkmatrices.c \
	checkmatrices_singleprec.c \
	nmf_als.c \
	nmf_als_singleprec.c \
	nmf_mu.c \
	nmf_mu_singleprec.c \
	outputtiming.c \
	setdefaultopts.c \
	setdefaultopts_singleprec.c \
        calculatenorm.c \
	calculatenorm_singleprec.c \
	generatematrix.c \
	generatematrix_singleprec.c \
	nmf_alspg.c \
	nmf_alspg_singleprec.c \
	nmf_neals.c \
	nmf_neals_singleprec.c \
	pg_subprob_h.c \
	pg_subprob_h_singleprec.c \
	pg_subprob_w.c \
	pg_subprob_w_singleprec.c \
	storematrix.c \
	storematrix_singleprec.c \
        calculatesvd.c \
	calculatesvd_singleprec.c
