/**
 * nmfDriver - performs a whole factorisation process including necessary steps before and after the factorisation
 * 
 * Purpose:
 *		This routine calculates a non negative matrix factorisation of a m by n matrix A
 *
 *		A = W * H
 *
 *		where A, W and H are non-negative matrices.
 *		W is a m by k matrix
 *		H is a k by n matrix
 *		k is the approximation factor
 *		
 * Description:
 *		This routine loads a matrix A (and in case they are specified matrices w0 and h0) from a file and
 *		performs multiple iterations of one of the implemented nmf-algorithms, until an iteration or tolerance
 *		limit is reached.
 *		The initial factor matrices W0 and H0 can be loaded from a file or initialised randomly.
 *		It is possible to run several complete runs, each with different inital matrices w0 and h0. For every run
 *		after the first the matrices then get randomly initialised.
 *		
 *		The factorisation results are normalised and rows of H and columns of W are resorted. The order is
 *		determined by the row sums of elementwise squared matrix H.
 *		
 *		Factorisation results of different runs are compared by the root mean square residual and the best result
 *		of all runs is stored in a file.
 *		Caution: The root mean square residual in this routine is the frobenius norm divided by sqrt(m*n)
 *
 *		Implemented algorithms to actually compute the factorisation are (for details see each algorithm):
 *
 *		nmf_mu		multiplicative update approach
 *		nmf_als		alternating least squares approach
 *		nmf_neals	normal equation alternating least squares approach
 *		nmf_alspg	alternating least squares approach using a projected gradient method
 *		nmf_pg		direct projected gradient approach
 *  * 		nmf_fast	Fast NMF algorithm
 * 		nmf_bayes	Bayes NMF algorithm
 *
 *
 *
 * Usage:
 *		For calling nmfDriver at least arguments a, k and iter have to be passed.
 *		Passing a NULL-pointer to w0 or h0 will let w0 and h0 be initialised according to opts->init
 *		Passing a NULL-pointer to opts will set default options
 *
 *		The datatype "options_t" has following structure:
 *
 * rep		in, 	number of repeated factorisations with differently initialised matrices w0 and h0
 *
 *
 * init		in, 	defines how to initialise the matrices w0 and h0
 *
 * min_init	in, 	defines the minimum value for initialisation
 *
 * max_init	in, 	defines the maximum value for initialisation
 *
 * w_out	in, 	filename to store final matrix w to
 *
 * h_out	in, 	filename to store final matrix h to
 *
 * TolX		in, 	tolerance value for convergence check of maxChange
 *			Every iteration checks for the change in each factor matrix
 * TolFun	in, 	tolerance value for convergence check of root mean square residual
 *
 *
 * Arguments:
 *
 * a		in, 	filename to load matrix a from
 *
 * k		in, 	approximation factor which limits the matrices dimensions
 *
 * iter		in, 	maximal number of iterations to perform
 *
 * w0		in, 	filename to load matrix w0 from, when NULL matrix is initialised using method in "opts->init"
 *
 * h0		in, 	filename to load matrix h0 from, when NULL matrix is initialised using method in "opts->init"
 *
 * alg		in,	algorithm to use for the factorisation steps
 *
 * opts		in,	options_t structure, which defines all additional options
 *			whenn NULL, defaultvalues are used
 */



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <sys/time.h>
#include <time.h>


#include "common.h"
#include "blaslapack.h"
#include "nmf_mu_singleprec.h"
#include "nmf_als_singleprec.h"
#include "nmf_neals_singleprec.h"
#include "nmf_alspg_singleprec.h"
#include "nmf_pg_singleprec.h"
#include "nmf_bayes_singleprec.h"
#include "nmf_fast_singleprec.h"
#include "checkmatrices_singleprec.h"
#include "loadmatrix_singleprec.h"
#include "storematrix_singleprec.h"
#include "generatematrix_singleprec.h"
#include "outputtiming.h"
#include "checkarguments_singleprec.h"
#include "setdefaultopts_singleprec.h"
#include "calculatenorm_singleprec.h"
#include "cpmidx_float.h"



void nmfDriver_singleprec(const char* a, const int k, int iter, const char* w0, const char* h0, alg_t alg, options_t_singleprec * opts) {


#ifdef PROFILE_NMF_DRIVER
	struct timeval start, end;
	gettimeofday(&start, 0);
#endif
#if DEBUG_LEVEL >= 2
	printf("Entering nmfDriver\n");
#endif

#ifdef ERROR_CHECKING
  errno = 0;		//no error occured so far
#endif

  options_t_singleprec default_opts;

  //Explicit options or Defaultoptions (if opts = NULL)
  if (!opts) {
    set_default_opts_singleprec(&default_opts);
    opts = &default_opts;
  }

  // checking arguments
  if (checkArguments_singleprec(a, k, iter, w0, h0, opts)) {
    perror("Arguments invalid in nmfDriver");
    return;
  }
  

  // declaring basic variables to load matrices
  int m, n, m_w0, n_w0, m_h0, n_h0;
  float * A, *W0, *H0;
  W0 = NULL;
  H0 = NULL;

  loadMatrix_singleprec(a, &m, &n, &A);
  if (errno) {
    perror("Loading matrix A failed in nmfDriver");
    free(A);
    return;
  }
  if (k > m || k > n) {
    errno = EDOM;
    perror("Approximation factor k bigger than original matrix dimensions.");
    free(A);
    return;
  }
  if (w0 && h0) {
    loadMatrix_singleprec(w0, &m_w0, &n_w0, &W0);
    loadMatrix_singleprec(h0, &m_h0, &n_h0, &H0);
  }  
  else {
    generateMatrix_singleprec(m, n, k, opts->init, opts->min_init, opts->max_init, &W0, &H0, A, opts);
  }  

#ifdef ERROR_CHECKING
  if(errno) {
    perror("Loading matrix W0 failed in nmfDriver");
    perror("Loading matrix H0 failed in nmfDriver");
    free(H0);
    free(A);
    free(W0);
    return;
  }
#endif


#ifdef PROFILE_MATLAB_COMPARISON
	struct timeval start_matlab, end_matlab;
	gettimeofday(&start_matlab, 0);
#endif


  if (checkMatrices_singleprec(A, W0, H0, m, n, k)) {
    errno = EDOM;
    perror("Matrices not compatible in nmfDriver");
    free(A);
    free(W0);
    free(H0);
    return;
  }
  
  
  //memory for saving matrices w and h with the best norm
  float * wbest = (float*) malloc(sizeof(float)*m*k);
  float * hbest = (float*) malloc(sizeof(float)*k*n);
  //memory for normalizing final matrices wbest and hbest
  float * hlen = (float*) malloc(sizeof(float)*k);
  idx_float * wlen = (idx_float*) malloc(sizeof(idx_float)*k);	//stores the indices of column-sums as well for resorting wbest/hbest
  

#ifdef ERROR_CHECKING
  if(errno) {
    perror("Failed allocating memory in nmfDriver");
    free(A);
    free(W0);
    free(H0);
    free(wbest);
    free(hbest);
    free(hlen);
    free(wlen);
    return;
  }
#endif

  float norm = 0.0;
  norm = HUGE_VAL;
  float normbest = 0.0;
  normbest = HUGE_VAL;

  int repetitions;

  for (repetitions = 1; repetitions <= opts->rep; ++repetitions) {
    switch (alg) {
      case mu:	norm = nmf_mu_singleprec(A, &W0, &H0, m, n, k, &iter, opts->TolX, opts->TolFun);
		break;
      case als:	norm = nmf_als_singleprec(&A, &W0, &H0, m, n, k, &iter, opts->TolX, opts->TolFun);
		break;
      case neals: norm = nmf_neals_singleprec(&A, &W0, &H0, m, n, k, &iter, opts->TolX, opts->TolFun);
		break;
      case alspg:	norm = nmf_alspg_singleprec(&A, &W0, &H0, m, n, k, &iter, opts->TolX);
		break;
      case pg:	norm = nmf_pg_singleprec(&A, &W0, &H0, m, n, k, &iter, opts->TolX);
		break;
      case bayes: norm = nmf_bayes_singleprec(&A, &W0, &H0, m, n, k, &iter, opts->TolX, opts->TolFun);
		break;
      case fast: norm = nmf_fast_singleprec(&A, &W0, &H0, m, n, k, &iter, opts->TolX, opts->TolFun);
		break;		
    }
    //storing best matrices w and h so far, if the norm is the best so far
    if (norm < normbest) {
      normbest = norm;
      swap_singleprec(&wbest, &W0);
      swap_singleprec(&hbest, &H0);
    }
    //re-initialise the starting matrices w0 and h0, if there are more repetitions to be run
    if (repetitions < opts->rep) {
      generateMatrix_singleprec(m, n, k, init_rand, opts->min_init, opts->max_init, &W0, &H0, A, opts);
    }
  } //end of loop from 1 to rep


#if DEBUG_LEVEL >= 0
  //Output final norm
  printf("Final Norm: %.16f\n", normbest);
#endif
 

  //normalizing results
  //-------------------
  float temp;
  temp = 0.;
  int i, j;

  //calculating hlen
  for (i = 0; i<k; ++i) {
    temp = 0.;
    for (j = 0; j<n; ++j) {
      temp += pow(hbest[i + j*k], 2);
    }
    temp = sqrt(temp);
    
    if (temp == 0.) {
      hlen[i] = 1.;
      fprintf(stderr, "Warning: Matrix H doesn't have full rank\n");
    }
    else
      hlen[i] = temp;
    
  }
  
  //wbest = wbest .* hlen'
  for (i=0; i<m; ++i) {
    for (j=0; j<k; ++j) {
      wbest[i + j*m] *= hlen[j];
    }
  }
  //hbest = hbest ./ hlen

  for (j=0; j<n; ++j)
    for (i=0; i<k; ++i)
      hbest[i + j*k] /= hlen[i];

  //Calculating wlen for sorting columns of w and rows of h
  for(j = 0; j<k; ++j) {
    temp = 0;
    for(i=0; i<m; ++i) {
      temp += wbest[i+j*m] * wbest[i+j*m];
    }
    wlen[j].val = temp;
    wlen[j].idx = j;
  }
  
  //sort wlen in descending order
  int elementsize = 0;
  elementsize = sizeof(idx_float);
  qsort(wlen, k, elementsize, cmpidx_float);
  
  
  //resorting columns of wbest according to order in wlen
  for(j=0; j<k; ++j) {
    for(i=0; i<m; ++i) {
      W0[i + j*m] = wbest[i + wlen[j].idx*m];
    }
  }
  
  //resorting rows of hbest according to order in wlen
  for(j=0; j<n; ++j) {
    for(i=0; i<k; ++i) {
      H0[i + j*k] = hbest[wlen[i].idx + j*k]; 
    }
  }



#ifdef PROFILE_MATLAB_COMPARISON
	gettimeofday(&end_matlab, 0);
	outputTiming("Timing of MATLAB_COMPARISON:", start_matlab, end_matlab);
#endif  
  

  //storing final results in files
  if (opts->w_out)
    storeMatrix_singleprec(opts->w_out, m, k, W0);
  if (opts->h_out)
    storeMatrix_singleprec(opts->h_out, k, n, H0);


  //freeing memory of dynamic variables
  free(A);
  free(W0);
  free(H0);
  free(wbest);
  free(hbest);
  free(hlen);
  free(wlen);

#if DEBUG_LEVEL >= 2
	printf("Exiting nmfDriver\n");
#endif
#ifdef PROFILE_NMF_DRIVER
	gettimeofday(&end, 0);
	outputTiming("Timing:", start, end);
#endif

}
//end of nmfDriver
//----------------
