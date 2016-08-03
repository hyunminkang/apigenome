//blaslapack.h
//version 0.95
//------------

//Purpose:	



#ifndef UNIVIE_NMFLIB_BLASLAPACK_H
#define UNIVIE_NMFLIB_BLASLAPACK_H


#ifdef __cplusplus
extern "C" {
#endif


//declaration of the BLAS/LAPACK routines to be used
//Fortran functions are wrapped in a new c-function
//--------------------------------------------------

//dsytrf - lapack routines
//-------------------------------
//used to compute a bunc-kaufman factorization of a Matrix a
static inline double dsytrf(char uplo, int n, double * a, int lda, int * ipiv, double * work, int lwork, int * info) {
  extern double dsytrf_(char * uplo, int * n, double * a, int * lda, int * ipiv, double * work, int * lwork, int * info);
  return dsytrf_(&uplo, &n, a, &lda, ipiv, work, &lwork, info);
}



//dsytrs - lapack routines
//-------------------------------
//used to solve a linear system of eq. using a bunc-kaufman factorization of a Matrix a
static inline double dsytrs(char uplo, int n, int nrhs, double * a, int lda, int * ipiv, double * b, int ldb, int * info) {
  extern double dsytrs_(char * uplo, int * n, int * nrhs, double * a, int * lda, int * ipiv, double * b, int * ldb, int * info);
  return dsytrs_(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, info);
}


//dnscrs - sparskit routine
//---------------------------------
// used to convert a dense matrix to CSR format
static inline double dnscsr(int nrow, int ncol, int nzmax, double * dns, int ndns, double * a, int * ja, int * ia, int * ierr) {
  extern double dnscsr_(int * nrow, int * ncol, int * nzmax, double * dns, int * ndns, double * a, int * ja, int * ia, int * ierr);
  return dnscsr_(&nrow, &ncol, &nzmax, dns, &ndns, a, ja, ia, ierr);
}

//aplsb - sparskit routine
//---------------------------------
// used to add two sparse matrices
static inline double aplsb(int nrow, int ncol, double * a, int * ja, int * ia, double s, double * b, int * jb, int * ib, double * c, int * jc, int * ic, int nzmax, int * ierr) {
  extern double aplsb_(int * nrow, int * ncol, double * a, int * ja, int * ia, double *s, double * b, int * jb, int * ib, double * c, int * jc, int * ic, int * nzmax, int * ierr);
  return aplsb_(&nrow, &ncol, a, ja, ia, &s, b, jb, ib, c, jc, ic, &nzmax, ierr);
}

//csrdns - sparskit routine
//---------------------------------
// used to convert a matrix in CSR format to a dense matrix
static inline double csrdns(int nrow, int ncol, double * a, int * ja, int * ia, double * dns, int ndns, int * ierr) {
  extern double csrdns_(int * nrow, int * ncol, double * a, int * ja, int * ia, double * dns, int * ndns, int * ierr);
  return csrdns_(&nrow, &ncol,  a, ja, ia, dns, &ndns, ierr);
}


//filter - sparskit routine
//----------------------------------
// used to filter entries according to their magnitude (csr-format)
static inline double filter(int n, int job, double drptol, double * a, int * ja, int * ia, double * b, int * jb, int *ib, int len, int * ierr) {
  extern double filter_(int * n, int * job, double *drptol, double * a, int * ja, int * ia, double * b, int * jb, int * ib, int * len, int * ierr);
  return filter_(&n, &job, &drptol, a, ja, ia, b, jb, ib, &len, ierr);
}


//getelm - sparskit routine
//---------------------------------
// used to get element (i,j) from sparse matrix
static inline double getelm (int i,int j,double * a, int * ja,int * ia,int iadd,int sorted) {
  extern double getelm_(int * i, int * j, double * a, int * ja, int * ia, int * iadd, int * sorted);
  return getelm_(&i, &j, a, ja, ia, &iadd, &sorted);
}


//amask - sparskit routine
//----------------------------------
// used to extract only those elements of a csr matrix that are indexed a by a mask pattern (pairs of row and column index)
static inline double amask(int nrow, int ncol, double * a, int * ja, int * ia, int * jmask, int * imask, double * c, int * jc, int * ic, int * iw, int nzmax, int * ierr) {
  extern double amask_(int * nrow, int * ncol, double * a, int * ja, int * ia, int * jmask, int * imask, double * c, int * jc, int * ic, int * iw, int * nzmax, int * ierr);
  return amask_(&nrow, &ncol, a, ja, ia, jmask, imask, c, jc, ic, iw, &nzmax, ierr);
}

//copmat - sparskit routine
//----------------------------
// used to copy a csr matrix to another csr  matrix
static inline double copmat(int nrow, double * a, int * ja, int * ia, double * ao, int * jao, int * iao, int ipos, int job) {
  extern double copmat_(int * nrow, double * a, int * ja, int * ia, double * ao, int *jao, int * iao, int * ipos, int * job);
  return copmat_(&nrow, a, ja, ia, ao, jao, iao, &ipos, &job);
}


//amub - sparskit routine
//-------------------------
// used to calculate C = A * B (all matrices in csr format)
static inline double amub(int nrow, int ncol, int job, double * a, int * ja, int * ia, double * b, int * jb, int * ib, double * c, int * jc, int * ic, int nzmax, int * iw, int * ierr) {
  extern double amub_(int * nrow, int * ncol, int * job, double * a, int * ja, int * ia, double * b, int * jb, int * ib, double * c, int * jc, int * ic, int * nzmax, int * iw, int * ierr);
  return amub_(&nrow, &ncol, &job, a, ja, ia, b, jb, ib, c, jc, ic, &nzmax, iw, ierr);
}
 


//transp - sparskit routine
//-----------------------------
// used to transpose a CSR matrix in place
static inline double transp(int nrow, int * ncol, double * a, int * ja, int * ia, int * iwk, int * ierr) {
  extern double transp_(int * nrow, int * ncol, double * a, int * ja, int * ia, int * iwk, int * ierr);
  return transp_(&nrow, ncol, a, ja, ia, iwk, ierr);
}


//csrcsc - sparskit routine
//------------------------------
// used to transpose a CSR matrix into CSC format (not in place)
static inline double csrcsc(int n, int job, int ipos, double * a, int * ja, int * ia, double * ao, int *jao, int * iao) {
  extern double csrcsc_(int * n, int * job, int * ipos, double * a, int * ja, int * ia, double * ao, int *jao, int * iao);
  return csrcsc_(&n, &job, &ipos, a, ja, ia, ao, jao, iao);
}
 
 
static inline double csrcsc2(int n, int n2, int job, int ipos, double * a, int * ja, int * ia, double * ao, int *jao, int * iao) {
  extern double csrcsc2_(int * n, int * n2, int * job, int * ipos, double * a, int * ja, int * ia, double * ao, int *jao, int * iao);
  return csrcsc2_(&n, &n2, &job, &ipos, a, ja, ia, ao, jao, iao);
}

//cg - sparskit routine
//-----------------------------
// iterative conjugate gradient method
static inline double cg(int n, double * rhs, double * sol, int * ipar, double * fpar, double * w) {
  extern double cg_(int * n, double * rhs, double * sol, int * ipar, double * fpar, double * w);
  return cg_(&n, rhs, sol, ipar, fpar, w);
}


//bcgstab - sparskit routine
//-----------------------------
// iterative conjugate gradient method
static inline double bcgstab(int n, double * rhs, double * sol, int * ipar, double * fpar, double * w) {
  extern double bcgstab_(int * n, double * rhs, double * sol, int * ipar, double * fpar, double * w);
  return bcgstab_(&n, rhs, sol, ipar, fpar, w);
}

//gmres - sparskit routine
//-----------------------------
// iterative generalized minimum residual method
static inline double gmres(int n, double * rhs, double * sol, int * ipar, double * fpar, double * w) {
  extern double gmres_(int * n, double * rhs, double * sol, int * ipar, double * fpar, double * w);
  return gmres_(&n, rhs, sol, ipar, fpar, w);
}

//cg - sparskit routine
//-----------------------------
// iterative conjugate gradient method
static inline double bcg(int n, double * rhs, double * sol, int * ipar, double * fpar, double * w) {
  extern double bcg_(int * n, double * rhs, double * sol, int * ipar, double * fpar, double * w);
  return bcg_(&n, rhs, sol, ipar, fpar, w);
}


//tfqmr - sparskit routine
//-----------------------------
// iterative conjugate gradient method
static inline double tfqmr(int n, double * rhs, double * sol, int * ipar, double * fpar, double * w) {
  extern double tfqmr_(int * n, double * rhs, double * sol, int * ipar, double * fpar, double * w);
  return tfqmr_(&n, rhs, sol, ipar, fpar, w);
}


//amux - sparskit routine
//--------------------------------
// sparse matrix * vector product
static inline double amux(int n, double * x, double * y, double * a, int * ja, int * ia) {
  extern double amux_(int * n, double * x, double * y, double * a, int * ja, int * ia);
  return amux_(&n, x, y, a, ja, ia);
}


//atmux - sparskit routine
//--------------------------------
// sparse transposed matrix * vector product
static inline double atmux(int n, double * x, double * y, double * a, int * ja, int * ia) {
  extern double atmux_(int * n, double * x, double * y, double * a, int * ja, int * ia);
  return atmux_(&n, x, y, a, ja, ia);
}

//ilut - sparskit routine
//----------------------------------
// incomplete LU factorization with dual truncation strategy -- preconditioner
static inline double ilut(int n, double * a, int * ja, int * ia, int lfil, double droptol, double * alu, int * jlu, int * ju, int iwk, double * w, int * jw, int * ierr) {
  extern double ilut_(int * n, double * a, int *  ja, int * ia, int *lfil, double * droptol, double * alu, int * jlu, int * ju, int * iwk, double * w, int * jw, int * ierr);
  return ilut_(&n, a, ja, ia, &lfil, &droptol, alu, jlu, ju, &iwk, w, jw, ierr);
}


static inline void ilutp(int n, double * a, int * ja, int * ia, int lfil, double droptol, double permtol, int mbloc, double * alu, int * jlu, int * ju, int iwk, double * w, int * jw, int * iperm, int *ierr) {
  extern void ilutp_(int * n, double * a, int * ja, int * ia, int * lfil, double * droptol, double *permtol, int * mbloc, double * alu, int * jlu, int * ju, int  *iwk, double *w, int * jw, int * iperm, int  *ierr);
  return ilutp_(&n, a, ja, ia, &lfil, &droptol, &permtol, &mbloc, alu, jlu, ju, &iwk, w, jw, iperm, ierr);
}


static inline void ilud(int n, double * a, int * ja, int * ia, double alpha, double droptol, double * alu, int * jlu, int * ju, int iwk, double * w, int * jw, int * ierr) {
  extern void ilud_(int * n, double * a, int * ja, int * ia, double *alpha, double *droptol, double * alu, int *jlu, int * ju, int * iwk, double * w, int * jw, int * ierr);  
  return ilud_(&n, a, ja, ia, &alpha, &droptol, alu, jlu, ju, &iwk, w, jw, ierr);
}


static inline void iludp(int n, double * a, int * ja, int * ia, double alpha, double droptol, double permtol, int mbloc, double * alu, int * jlu, int * ju, int iwk, double * w, int * jw, int * iperm, int * ierr) {
  extern void iludp_(int * n, double * a, int * ja, int * ia, double *alpha, double *droptol, double  * permtol, int *mbloc, double * alu, int *jlu, int * ju, int * iwk, double * w, int * jw, int * iperm, int * ierr);  
  return iludp_(&n, a, ja, ia, &alpha, &droptol, &permtol, &mbloc, alu, jlu, ju, &iwk, w, jw, iperm, ierr);
}



//lusol -- sparkit routine
//------------------------------
// solves the system (LU) x = y for an LU decomposition stored in (alu, jlu, ju) msr-format
static inline double lusol(int n, double * y, double * x, double * alu, int * jlu, int * ju) {
  extern double lusol_(int * n, double * y, double * x, double * alu, int * jlu, int * ju);
  return lusol_(&n, y, x, alu, jlu, ju);
}

//lutsol -- sparkit routine
//-------------------------------
// solves the system Transp(LU) x = y for an LU decomposition stored in (alu, jlu, ju) msr-format
static inline double lutsol(int n, double * y, double * x, double * alu, int * jlu, int * ju) {
  extern double lutsol_(int * n, double * y, double * x, double * alu, int * jlu, int * ju);
  return lutsol_(&n, y, x, alu, jlu, ju);
}



// dtrsm - Blas level routine
//----------------------------------
// used to solve a linear system A * x = B with A beeing triangular
static inline double dtrsm(char side, char uplo, char transa, char diag, int m, int n, double alpha, double * a, int lda, double * b, int ldb) {
  extern double dtrsm_(char* side, char* uplo, char* transa, char* diag, int* m, int* n, double* alpha, double * a, int* lda, double * b, int* ldb);
  return dtrsm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
}

static inline float strsm(char side, char uplo, char transa, char diag, int m, int n, float alpha, float * a, int lda, float * b, int ldb) {
extern float strsm_(char* side, char* uplo, char* transa, char* diag, int* m, int* n, float* alpha, float * a, int* lda, float * b, int* ldb);
return strsm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
}

// dlange - Lapack auxiliary routine
//----------------------------------
// used to calculate the absolute maximum of a matrix and to calculate the frobenius norm
static inline double dlange(char norm, int m, int n, double* a, int lda, double * work) {
  extern double dlange_(char* norm, int* m, int* n, double* a, int* lda, double* work);
  return dlange_(&norm, &m, &n, a, &lda, work);
}

// slange - LAPACK auxiliary routine
//-----------------------------------
// used to calculate the absolut maximum of a matrix and to calculate the frobenius norm
static inline float slange(char norm, int m, int n, float* a, int lda, float * work) {
  extern float slange_(char* norm, int* m, int* n, float* a, int* lda, float* work);
  return slange_(&norm, &m, &n, a, &lda, work);
}


// dlamch - Lapack auxiliary routine
//----------------------------------
// used to calculate machine precision epsilon
static inline double dlamch(char cmach) {
  extern double dlamch_(char* cmach);
  return dlamch_(&cmach);
}

// slamch - Lapack auxiliary routine - single precision !!!
//---------------------------------------------------------
// used to calculate machine precision epsilon
static inline float slamch(char cmach) {
  extern float slamch_(char* cmach);
  return slamch_(&cmach);
}


// dgemm - BLAS routine
//---------------------
// used to calculate general matrix-matrix multiplications of the form
// C = alpha * A * B + beta * C
static inline int dgemm(char transa, char transb, int m, int n, int k, double alpha, double* a, int lda, double* b, int ldb, double beta, double* c, int ldc) {
  extern int dgemm_(char* transa, char* transb, int* m, int* n, int* k, double* alpha, double* a, int* lda, double* b, int* ldb, double* beta, double* c, int* ldc);
  return dgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

// sgemm - BLAS routine - single precision !!!
//---------------------------------------------
// used to calculate general matrix-matrix multiplications of the form
// C = alpha * A * B + beta * C
static inline int sgemm(char transa, char transb, int m, int n, int k, float alpha, float* a,  int lda, float* b, int ldb, float beta, float* c, int ldc) {
  extern int sgemm_(char* transa, char* transb, int* m, int* n, int* k, float* alpha, float* a, int* lda, float* b, int* ldb, float* beta, float* c, int* ldc);
  return sgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}



// dgeqp3 - Lapack routine
//------------------------
// used to calculate a QR factorization with column pivoting
static inline int dgeqp3(int M, int N, double *A, int LDA, int * JPVT, double * TAU, double * WORK, int LWORK, int *INFO ) {
  extern int dgeqp3_(int *M, int *N, double *A, int *LDA, int *JPVT, double *TAU, double *WORK, int * LWORK, int *INFO);
  return dgeqp3_(&M, &N, A, &LDA, JPVT, TAU, WORK, &LWORK, INFO);
}

static inline int sgeqp3(int M, int N, float *A, int LDA, int * JPVT, float * TAU, float * WORK, int LWORK, int *INFO ) {
  extern int sgeqp3_(int *M, int *N, float *A, int *LDA, int *JPVT, float *TAU, float *WORK, int * LWORK, int *INFO);
  return sgeqp3_(&M, &N, A, &LDA, JPVT, TAU, WORK, &LWORK, INFO);
}


// dormqr - Lapack routine
//-----------------------
// used to calculate a matrix-matrix-multiplication of a factor matrix q (result of dgeqp3) and another matrix
static inline int dormqr(char SIDE, char TRANS, int M, int N, int K, double * A, int LDA, double * TAU, double * C, int LDC, double * WORK, int LWORK, int *INFO ) {
  extern int dormqr_(char * SIDE, char * TRANS, int * M, int * N, int * K, double * A, int * LDA, double * TAU, double * C, int * LDC, double * WORK, int * LWORK, int * INFO);
  return dormqr_(&SIDE, &TRANS, &M, &N, &K, A, &LDA, TAU, C, &LDC, WORK, &LWORK, INFO);
}


static inline int sormqr(char SIDE, char TRANS, int M, int N, int K, float * A, int LDA, float * TAU, float * C, int LDC, float * WORK, int LWORK, int *INFO ) {
  extern int sormqr_(char * SIDE, char * TRANS, int * M, int * N, int * K, float * A, int * LDA, float * TAU, float * C, int * LDC, float * WORK, int * LWORK, int * INFO);
  return sormqr_(&SIDE, &TRANS, &M, &N, &K, A, &LDA, TAU, C, &LDC, WORK, &LWORK, INFO);
}


// dtrtrs - Lapack routine
//------------------------
// used to solve a triangular system of the form a * x = b
static inline int dtrtrs(char UPLO, char TRANS, char DIAG, int N, int NRHS, double * A, int LDA, double * B, int LDB, int *INFO) {
  extern int dtrtrs_(char *UPLO, char *TRANS, char *DIAG, int *N, int *NRHS, double * A, int *LDA, double * B, int *LDB, int *INFO);
  return dtrtrs_(&UPLO, &TRANS, &DIAG, &N, &NRHS, A, &LDA, B, &LDB, INFO);
}


static inline int strtrs(char UPLO, char TRANS, char DIAG, int N, int NRHS, float * A, int LDA, float * B, int LDB, int *INFO) {
  extern int strtrs_(char *UPLO, char *TRANS, char *DIAG, int *N, int *NRHS, float * A, int *LDA, float * B, int *LDB, int *INFO);
  return strtrs_(&UPLO, &TRANS, &DIAG, &N, &NRHS, A, &LDA, B, &LDB, INFO);
}


// dgesv - Lapack routine
//-----------------------
// used to calculate a LU-Factorisation of a matrix
static inline int dgesv(int N, int NRHS, double * A, int LDA, int * IPIV, double * B, int LDB, int * INFO) {
  extern int dgesv_(int * N, int * NRHS, double * A, int * LDA, int * IPIV, double * B, int * LDB, int * INFO);
  return dgesv_(&N, &NRHS, A, &LDA, IPIV, B, &LDB, INFO);
}


static inline int sgesv(int N, int NRHS, float * A, int LDA, int * IPIV, float * B, int LDB, int * INFO) {
  extern int sgesv_(int * N, int * NRHS, float * A, int * LDA, int * IPIV, float * B, int * LDB, int * INFO);
  return sgesv_(&N, &NRHS, A, &LDA, IPIV, B, &LDB, INFO);
}

// dlacpy - Lapack routine
//-------------------------
// used to copy a matrix A or part of it to a matrix B
static inline int dlacpy(char UPLO, int M, int N, double * A, int LDA, double * B, int LDB) {
  extern int dlacpy_(char * UPLO, int * M, int * N, double * A, int * LDA, double * B, int * LDB);
  return dlacpy_(&UPLO, &M, &N, A, &LDA, B, &LDB);
}

//slacpy - LAPACK routine
//-----------------------
// used to copy a matrix A or part of it to a matrix B
static inline int slacpy(char UPLO, int M, int N, float * A, int LDA, float * B, int LDB) {
  extern int slacpy_(char * UPLO, int* M, int* N, float* A, int* LDA, float* B, int* LDB);
  return slacpy_(&UPLO, &M, &N, A, &LDA, B, &LDB);
}



// dscal - Blas routine
//---------------------
// used for multiplying a vector by a constant
static inline int dscal(int N, double DA, double * DX, int INCX) {
  extern int dscal_(int *N, double * DA, double * DX, int * INCX);
  return dscal_(&N, &DA, DX, &INCX);
}

static inline int sscal(int N, float DA, float * DX, int INCX) {
  extern int sscal_(int *N, float * DA, float * DX, int * INCX);
  return sscal_(&N, &DA, DX, &INCX);
}


// dcopy - Blas routine
//---------------------
// used for copying a vector
static inline int dcopy(int N, double * DX, int INCX, double * DY, int INCY) {
  extern int dcopy_(int * N, double * DX, int * INCX, double *DY, int * INCY);
  return dcopy_(&N, DX, &INCX, DY, &INCY);
}

static inline int scopy(int N, float * DX, int INCX, float * DY, int INCY) {
  extern int scopy_(int * N, float * DX, int * INCX, float *DY, int * INCY);
  return scopy_(&N, DX, &INCX, DY, &INCY);
}

// daxpy - Blas level 1 routine
//-----------------------------
// used for calculating dy = da * dx + dy
static inline int daxpy(int N, double DA, double * DX, int INCX, double * DY, int INCY ) {
  extern int daxpy_(int * N, double * DA, double * DX, int * INCX, double * DY, int *INCY);
  return daxpy_(& N, & DA, DX, &INCX, DY, &INCY);
}


// saxpy - BLAS level 1 routine
//-----------------------------
// used for calculating dy = da * dx + dy
static inline int saxpy(int N, float DA, float* DX, int INCX, float* DY, int INCY) {
  extern int saxpy_(int * N, float * DA, float* DX, int* INCX, float* DY, int* INCY);
  return saxpy_(&N, &DA, DX, &INCX, DY, &INCY);  
}


// dorgqr - Lapack routine
//------------------------
// used for calculating the explicit matrix Q out of a matrix factorized and overwritten by dgeqp3
static inline int dorgqr(int M, int N, int K, double * A, int LDA, double * TAU, double * WORK, int LWORK, int * INFO) {
  extern int dorgqr_(int * M, int * N, int * K, double * A, int * LDA, double * TAU, double * WORK, int * LWORK, int * INFO);
  return dorgqr_(&M, &N, &K, A, &LDA, TAU, WORK, &LWORK, INFO);
}



static inline int sorgqr(int M, int N, int K, float  * A, int LDA, float  * TAU, float  * WORK, int LWORK, int * INFO) {
  extern int sorgqr_(int * M, int * N, int * K, float  * A, int * LDA, float  * TAU, float  * WORK, int * LWORK, int * INFO);
  return sorgqr_(&M, &N, &K, A, &LDA, TAU, WORK, &LWORK, INFO);
}

// PROPACK SVD routine
static inline int dlansvd(char jobu, char jobv, int m, int n, int k, int kmax, int (*aprod)(char*, int*, int*, double*, double*, double*, int*), double* U, int ldu, double * Sigma, double* bnd, double* V, int ldv, double tolin, double * work, int lwork, int* iwork, int liwork, double * doption, int * ioption, int * info, double * dparm, int * iparm) {
    extern int dlansvd_(char * jobu, char * jobv, int* m, int* n, int * k, int* kmax, int (*aprod)(char*, int*, int*, double*, double*, double*, int*),double * U, int * ldu, double * Sigma, double * bnd, double * V, int * ldv,double * tolin, double * work, int * lwork, int *iwork, int * liwork,double * doption,int * ioption,int * info,double * dparm, int * iparm);
    return dlansvd_(&jobu, &jobv, &m, &n, &k, &kmax, aprod, U, &ldu, Sigma, bnd, V, &ldv, &tolin, work, &lwork, iwork, &liwork, doption, ioption, info, dparm, iparm); 
}

static inline int slansvd(char jobu, char jobv, int m, int n, int k, int kmax, int (*aprod)(char*, int*, int*, float *, float *, float *, int*), float * U, int ldu, float  * Sigma, float * bnd, float * V, int ldv, float  tolin, float  * work, int lwork, int* iwork, int liwork, float  * doption, int * ioption, int * info, float  * dparm, int * iparm) {
    extern int slansvd_(char * jobu, char * jobv, int* m, int* n, int * k, int* kmax, int (*aprod)(char*, int*, int*, float *, float *, float *, int*),float  * U, int * ldu, float  * Sigma, float  * bnd, float  * V, int * ldv,float  * tolin, float  * work, int * lwork, int *iwork, int * liwork,float  * doption,int * ioption,int * info,float  * dparm, int * iparm);
    return slansvd_(&jobu, &jobv, &m, &n, &k, &kmax, aprod, U, &ldu, Sigma, bnd, V, &ldv, &tolin, work, &lwork, iwork, &liwork, doption, ioption, info, dparm, iparm); 
}

// dgemv - Lapack routine
//-----------------------
// used to calculate one of the matrix-vector-operations y := alpha * A * x + beta * y _or_ y := alpha * A' * x + beta * y
static inline int dgemv(char TRANS, int M, int N, double ALPHA, double * A, int LDA, double * X, int INCX, double BETA, double * Y, int INCY) {
  extern int dgemv_(char * TRANS, int * M, int * N, double * ALPHA, double * A, int * LDA, double * X, int * INCX, double * BETA, double * Y, int * INCY);
  return dgemv_(&TRANS, &M, &N, &ALPHA, A, &LDA, X, &INCX, &BETA, Y, &INCY);
}

static inline int sgemv(char TRANS, int M, int N, float  ALPHA, float  * A, int LDA, float  * X, int INCX, float  BETA, float  * Y, int INCY) {
  extern int sgemv_(char * TRANS, int * M, int * N, float  * ALPHA, float  * A, int * LDA, float  * X, int * INCX, float  * BETA, float  * Y, int * INCY);
  return sgemv_(&TRANS, &M, &N, &ALPHA, A, &LDA, X, &INCX, &BETA, Y, &INCY);
}


// dlaset - Lapack routine
//------------------------
// used to initialize a m by n matrix to a fixed scalar
static inline int dlaset(char uplo, int m, int n, double alpha, double beta, double * A, int lda) {
  extern int dlaset_(char * uplo, int * m, int * n, double * alpha, double * beta, double * A, int * lda);
  return dlaset_(&uplo, &m, &n, &alpha, &beta, A, &lda);
}



static inline int slaset(char uplo, int m, int n, float  alpha, float  beta, float  * A, int lda) {
  extern int slaset_(char * uplo, int * m, int * n, float  * alpha, float  * beta, float  * A, int * lda);
  return slaset_(&uplo, &m, &n, &alpha, &beta, A, &lda);
}


// dsaupd  ARPACK reverse communication interface routine.
//--------------------------------------------------------
static inline int dsaupd(int * IDO, unsigned char BMAT, int N, unsigned char * WHICH, int NEV, double TOL, double * RESID, int NCV, double * V, int LDV, int * IPARAM, int * IPNTR, double * WORKD, double * WORKL, int LWORKL, int * INFO) {
  extern int dsaupd_(int * IDO, unsigned char * BMAT, int * N, unsigned char * WHICH, int * NEV, double * TOL, double * RESID, int * NCV, double * V, int * LDV, int * IPARAM, int * IPNTR, double * WORKD, double * WORKL, int * LWORKL, int * INFO);
  return dsaupd_(IDO, &BMAT, &N, WHICH, &NEV, &TOL, RESID, &NCV, V, &LDV, IPARAM, IPNTR, WORKD, WORKL, &LWORKL, INFO );
}


static inline int ssaupd(int * IDO, unsigned char BMAT, int N, unsigned char * WHICH, int NEV, float  TOL, float  * RESID, int NCV, float  * V, int LDV, int * IPARAM, int * IPNTR, float  * WORKD, float  * WORKL, int LWORKL, int * INFO) {
  extern int ssaupd_(int * IDO, unsigned char * BMAT, int * N, unsigned char * WHICH, int * NEV, float  * TOL, float  * RESID, int * NCV, float  * V, int * LDV, int * IPARAM, int * IPNTR, float  * WORKD, float  * WORKL, int * LWORKL, int * INFO);
  return ssaupd_(IDO, &BMAT, &N, WHICH, &NEV, &TOL, RESID, &NCV, V, &LDV, IPARAM, IPNTR, WORKD, WORKL, &LWORKL, INFO );
}

// dseupd  ARPACK routine that returns Ritz values and (optionally) Ritz vectors
//----------------------------------------------------------------
static inline int dseupd(int RVEC, char HOWMNY, int * SELECT, double * D, double * Z, int LDZ, double SIGMA, unsigned char BMAT, int N, unsigned char * WHICH, int NEV, double TOL, double * RESID, int NCV, double * V, int LDV, int * IPARAM, int * IPNTR, double * WORKD, double *  WORKL, int LWORKL, int * INFO ) {
  extern int dseupd_(int * RVEC, char * HOWMNY, int * SELECT, double * D, double * Z, int * LDZ, double * SIGMA, unsigned char *BMAT, int *N, unsigned char * WHICH, int * NEV, double * TOL, double * RESID, int * NCV, double * V, int *LDV, int * IPARAM, int * IPNTR, double * WORKD, double *  WORKL, int * LWORKL, int * INFO);
  return dseupd_(&RVEC, &HOWMNY, SELECT, D, Z, &LDZ, &SIGMA, &BMAT, &N, WHICH, &NEV, &TOL, RESID, &NCV, V, &LDV, IPARAM, IPNTR, WORKD, WORKL, &LWORKL, INFO);
}  

static inline int sseupd(int RVEC, char HOWMNY, int * SELECT, float  * D, float  * Z, int LDZ, float  SIGMA, unsigned char BMAT, int N, unsigned char * WHICH, int NEV, float  TOL, float  * RESID, int NCV, float  * V, int LDV, int * IPARAM, int * IPNTR, float  * WORKD, float  *  WORKL, int LWORKL, int * INFO ) {
  extern int sseupd_(int * RVEC, char * HOWMNY, int * SELECT, float  * D, float  * Z, int * LDZ, float  * SIGMA, unsigned char *BMAT, int *N, unsigned char * WHICH, int * NEV, float  * TOL, float  * RESID, int * NCV, float  * V, int *LDV, int * IPARAM, int * IPNTR, float  * WORKD, float  *  WORKL, int * LWORKL, int * INFO);
  return sseupd_(&RVEC, &HOWMNY, SELECT, D, Z, &LDZ, &SIGMA, &BMAT, &N, WHICH, &NEV, &TOL, RESID, &NCV, V, &LDV, IPARAM, IPNTR, WORKD, WORKL, &LWORKL, INFO);
}  


// dnrm2 - Blas routine
//---------------------
// used for calculating the euclidean norm of a vector x
static inline double dnrm2(int N, double * X, int INCX) {
  extern double dnrm2_(int *N, double * X, int * INCX);
  return dnrm2_(&N, X, &INCX);
}

static inline float  snrm2(int N, float  * X, int INCX) {
  extern float  snrm2_(int *N, float  * X, int * INCX);
  return snrm2_(&N, X, &INCX);
}

// dswap BLAS level 1 routine to swap to vectors
//-----------------------------------------------
static inline int dswap(int N, double * DX, int INCX, double * DY, int INCY) {
  extern int dswap_(int * N, double * DX, int * INCX, double * DY, int *INCY);
  return dswap_(&N, DX, &INCX, DY, &INCY);
} 


static inline int sswap(int N, float  * DX, int INCX, float  * DY, int INCY) {
  extern int sswap_(int * N, float  * DX, int * INCX, float  * DY, int *INCY);
  return sswap_(&N, DX, &INCX, DY, &INCY);
} 


// ddot -- BLAS routine to calculate the dot product -- ATTENTION -- called distdot_ because SPARSPACK solver need this routine with the same signature
static inline double distdot_(int N, double * DX, int INCX, double * DY, int INCY) {
	extern int ddot_(int * N, double * DX, int * INCX, double * DY, int * INCY);
	return ddot_(&N, DX, &INCX, DY, &INCY);
}

#ifdef __cplusplus
}
#endif


#endif

