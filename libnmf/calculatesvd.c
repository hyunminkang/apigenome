#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>

#include "common.h"
#include "outputtiming.h"
#include "blaslapack.h"


/**
 * Calculates the Singular value decomposition for a matrix a
 * 
 * Purpose	
 *		Calculates the SVD for a matrix a, needed for the NNDSVD initialization.
 *		@see generateMatrices
 * 
 * Description
 *		Uses ARPACK routines dsaupd/dseupd for calcuation.
 * 
 * Arguments
 *
 *	a			in,	pointer to matrix to be factorized a = u * diag(s) * v'
 *	u			out,	pointer to factor matrix u
 *	s			out,	pointer to vector of singular values
 *	v			out,	pointer to factor matrix v
 *	m			in, 	rows of a
 *	n			in, 	cols of a
 * 	k			in, 	rank of approximation
 *	maxiter			in,	maximum number of iterations
 *	tol			in,	tolerance limit
 *	ncv			in,	length of arnoldi factorization
 *	order			in,	descending order of singular values if order = 1, ascending if order = 0
 *
 *	return value		number of converged ritz values -- on exit compare with original k
 */
int calculateSVD(double * a, double * u, double * s, double * v, int m, int n, int k, int maxiter, double tol, int ncv, int order)
{
#ifdef PROFILE_CALCULATE_SVD
	struct timeval start, end;
	gettimeofday(&start, 0);
#endif
#if DEBUG_LEVEL >= 2
	printf("Entering calculateSVD\n");
#endif

  //Matrix-Dimensions
  int mindim, maxdim;
  mindim = (m < n) ? m : n;
  maxdim = (m > n) ? m : n;
  //Parameter tests
  if (m < 1 || n < 1)				//illegal matrix dimensions
    return -1;
  if (k < 1 || k >= mindim)		        //illegal approximation factor: 0 < k < mindim
    return -2;				
  if (maxiter < 1 || tol < 0)			//illegal tolerance or maxiter parameter
    return -3;
  if (ncv > mindim || ncv <= k)			//illegal ncv value: k < ncv <= mindim
    return -4;
    
    
  //Reset errno to 0 - no error yet
  errno = 0;
  

  
  int ldv;
  //local arrays
  double * workl, *workd, *resid, *ax;
  int * select;
  int * iparam, * ipntr;
  iparam = (int*) malloc(sizeof(int)* 11);
  ipntr = (int*) malloc(sizeof(int)*11);
  
  
  
  
  //local scalars
  unsigned char bmat;
  unsigned char which[2];
  which[0] = 'L';	//Ask for the NEV singular values of largest magnitude   
  which[1] = 'M';
  int ido, nev, info, ierr, nconv, rvec;
  ierr = 0;
  double sigma, temp;
  sigma = 0;
  
  //trace information
  int  logfil, ndigit, msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd;
  
  ndigit = -3;
  logfil = 6;
  msgets = 0;
  msaitr = 0;
  msapps = 0;
  msaupd = 1;
  msaup2 = 1;
  mseigt = 0;
  mseupd = 1;
  
  //nev -> number of singular values to be computed
  //ncv -> sets the length of the arnoldi factorization
  nev = k;
  ldv = n;
  
  
  //allocating memory for dynamic arrays

  workl = (double*) malloc(sizeof(double) * ncv * (ncv + 8));	//working space matrix
  workd = (double*) malloc(sizeof(double)*3*(mindim));			//working space matrix
  resid = (double*) malloc(sizeof(double)*maxdim);			//residual vector
  ax = (double*) malloc(sizeof(double)*maxdim);				// used for computing matrix vector products
  select = (int*) malloc(sizeof(int)*ncv);				// used for ritz vector permutation in dseupd
  
  

  bmat = 'I';		//standard problem B = I
  		
  
  
  // Specification of stopping rules initial conditions before calling DASUPD
  
  int lworkl = ncv * (ncv+8);
  // Initializing 
  info = 0;			//dsaupd status codes
  ido = 0;			//dsaupd user task to be done
  
  // Specification of algorithm mode
  //--------------------------------
  
  iparam[0] = 1;		//use exact shift strategy
  iparam[2] = maxiter;		//maximum number of Arnoldi iterations allowed
  iparam[6] = 1;		//use mode 1 of DSAUPD
  
  
  
  
  //main loop (reverse communication loop)

  while(1){

    //call dsaupd
    if (m >= n) {
      dsaupd(&ido, bmat, mindim, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl,lworkl, &info );
    }
    else {
      dsaupd(&ido, bmat, mindim, which, nev, tol, resid, ncv, u, m, iparam, ipntr, workd, workl,lworkl, &info );
    }
    
    // error or finished?
    if (ido != -1 && ido != 1)
      break;
    
    //perform matrix vector multiplications
    if (m >= n) {
      dgemv('N', m, n, 1.0, a, m, workd + ipntr[0] - 1, 1, 0.0, ax, 1); 
      dgemv('T', m, n, 1.0, a, m, ax, 1, 0.0, workd + ipntr[1] - 1, 1);
    }
    else {
      dgemv('T', m, n, 1.0, a, m, workd + ipntr[0] - 1, 1, 0.0, ax, 1); 
      dgemv('N', m, n, 1.0, a, m, ax, 1, 0.0, workd + ipntr[1] - 1, 1);
    }
  }
  
  //either we have convergence or there is an error
  if (info < 0) {
    return -5;					//Error with _saupd
  }
  else {
    //no fatal errors occured -> post-processing using DSEUPD
    //compute all ritz vectors
    rvec = 1;
    if (m >= n) {
      dseupd( rvec, 'A', select, s, v, ldv, sigma, bmat, mindim, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, &ierr);
    }
    else {
      dseupd( rvec, 'A', select, s, u, m, sigma, bmat, mindim, which, nev, tol, resid, ncv, u, m, iparam, ipntr, workd, workl, lworkl, &ierr);
    }
  }

//        %-----------------------------------------------%
//        | Singular values are returned in the first     |
//        | column of the two dimensional array S         |
//        | and the corresponding right singular vectors  | 
//        | are returned in the first NEV columns of the  |
//        | two dimensional array V (m >= n) or U (m < n) |
//        %-----------------------------------------------%

  if (ierr) {
    return -6;					//Error with _seupd
  }
  else {
    // get number of converged ritz values
    nconv = iparam[4];
    int j;
    //for every converged ritz value
    for(j = 0; j < nconv; ++j) {
      s[j] = sqrt(s[j]);			//calculate
      
// 	      	 %-----------------------------%
//               | Compute the left singular   |
//               | vectors from the formula    |
//               |                             |
//               |     u = Av/sigma    (m >= n)|
//               |     v = A'u/sigma   (m < n) |
//               | u should have norm 1 so     |
//               | divide by norm(Av)/norm(Au) |
//               %-----------------------------%

      if (m >= n) {
        dgemv('N', m, n, 1.0, a, m, v + n*j, 1, 0.0, ax, 1);	//ax = A*v
        dcopy(maxdim, ax, 1, u + j*m, 1);
        temp = 1.0/dnrm2(maxdim, u + j*m, 1);
        dscal(maxdim, temp, u + j*m, 1);
      }
      else {
	dgemv('T', m, n, 1.0, a, m, u + m*j, 1, 0.0, ax, 1);	//ax = A'*u
        dcopy(maxdim, ax, 1, v + j*n, 1);
        temp = 1.0/dnrm2(maxdim, v + j*n, 1);
        dscal(maxdim, temp, v + j*n, 1);
      }
      
    }

  }  
  
  // Reordering of results if order = 1 (ascending) since s. values and vectors are returned in descending order by default
  if (order && (k > 1) ) {
    //reordering  matrix u
    int index;
    double tmp;
    for (index = 0; index < (k / 2); ++ index)
      dswap(m, u + index * m, 1, u + (k - index - 1) * m, 1);
    
    //reordering vector s
    for (index = 0; index < (k / 2); ++ index) {
      tmp = s[index];
      s[index] = s[k - index - 1];
      s[k - index - 1] = tmp;
    }
    
    //reordering matrix v
    for (index = 0; index < (nev / 2); ++ index)
      dswap(n, v + index * n, 1, v + (k - index - 1) * n, 1);    
  }

  
  
#if DEBUG_LEVEL >= 2
	printf("Exiting calculateSVD\n");
#endif
#ifdef PROFILE_CALCULATE_SVD
	gettimeofday(&end, 0);
	outputTiming("SVD Timing:", start, end);
#endif
 
  free(workl);
  free(workd);
  free(resid);
  free(ax);
  free(select);
  free(iparam);
  free(ipntr);
  return nconv;			//returning number of singular values which reached the given tolerance
}

