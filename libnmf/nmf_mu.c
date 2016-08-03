/**
 * NMF_MU -	calculates the nmf using a multiplicative update method
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
 *		The multiplicative update method as described in Berry (2006) (see references [1])
 *		is based on following algorithm (in Matlab notation):
 *
 *		W = rand(m, k);
 *		H = rand(k, n);
 *		for iter = 1:maxiter
 *			H = H .* (W' * A) ./ (W' * W * H + 10E-09)
 *			W = W .* (A * H') ./ (W * H * H' + 10E-09)
 *		end
 *
 *		This routine uses the BLAS "dgemm" routine for matrix-matrix-multiplications
 *
 * Arguments:
 *
 * a		in, 	pointer to matrix to factorise
 *
 * w0		in, 	pointer to initialised factor-matrix w0
 *		out, 	pointer to final factor-matrix w
 *
 * h0		in, 	pointer to initialised factor-matrix h0
 *		out, 	pointer to final factor-matrix h
 *
 * m		in, 	first dimension of a and w/w0
 *
 * n		in, 	second dimension of a and h/h0
 *
 * k		in, 	second dimension of w/w0 and first dimension of h/h0
 *
 * maxiter	in, 	maximal number of iterations to run
 *		out, 	number of iterations actually run
 *
 * TolX		in, 	used in check for convergence, tolerance for maxchange of matrices w and h
 *
 * TolFun	in, 	used in check for convergence, tolerance for root mean square residual
 */




//defines a factor added to matrix elements when used as divisor to avoid division by zero
//----------------------------------------------------------------------------------------
#define DIV_BY_ZERO_AVOIDANCE 1E-09






#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <sys/time.h>


#include "common.h"
#include "blaslapack.h"
#include "outputtiming.h"
#include "calculatenorm.h"
#include "calculatemaxchange.h"







double nmf_mu(double * a, double ** w0, double ** h0, int m, int n, \
		      int k, int * maxiter, const double TolX, const double TolFun) 
{
  

#ifdef PROFILE_NMF_MU
        struct timeval start, end;
        gettimeofday(&start, 0);
#endif

#if DEBUG_LEVEL >= 2
	printf("Entering nmf_mu\n");
#endif


#ifdef ERROR_CHECKING
  errno = 0;
#endif

  // definition of necessary dynamic data structures
  //...for calculating matrix h
  double* numerh = (double*) malloc(sizeof(double) *k*n);
  double* work1 = (double*) malloc(sizeof(double)*k*k);					// used for calculation of h & w
  double* work2 = (double*) malloc(sizeof(double)*k*n);
  double* h = (double*) malloc(sizeof(double)*k*n);
  //----------------

  //...for calculating matrix w
  double* numerw = (double*) malloc(sizeof(double)*m*k);
  double* work2w = (double*) malloc(sizeof(double)*m*k);
  double* w = (double*) malloc(sizeof(double)*m*k);
  //-----------------

  //...for calculating the norm of A-W*H
  double* d = (double*) malloc(sizeof(double)*m*n);					//d = a - w*h
  double dnorm0 = 0.;
  double dnorm = 0.;
  const double eps = dlamch('E');					//machine precision epsilon
  const double sqrteps = sqrt(eps);					//squareroot of epsilon
  


  //-------------------

#ifdef ERROR_CHECKING
  if (errno) {
    perror("Error allocating memory in nmf_mu");
    free(numerh);
    free(numerw);
    free(work1);
    free(work2);
    free(h);
    free(w);
    free(work2w);
    free(d);
    return -1;
  }
#endif



//Is ZERO_THRESHOLD _not_ defined then use machine epsilon as default
#ifndef ZERO_THRESHOLD
  const double ZERO_THRESHOLD = eps;
#endif

  
  //Loop-Indices
  int iter, i;



  // factorisation step in a loop from 1 to maxiter
  for (iter = 1; iter <= *maxiter; ++iter) {



    // calculating matrix h
    //----------------
    // calculating numerh = w0'*a
    dgemm('T', 'N', k, n, m, 1.0, *w0, m, a, m, 0., numerh, k);  
    // calculating first intermediate result work1 = w0'*w0
    dgemm('T', 'N', k, k, m, 1.0, *w0, m, *w0, m, 0., work1, k);
    // calculating second intermediate result work2 = work1 * h0
    dgemm('N', 'N', k, n, k, 1.0, work1, k, *h0, k, 0., work2, k);


    //calculating elementwise matrixmultiplication, Division and addition h = h0 .* (numerh ./(work2 + eps))
    //set elements < zero_threshold to zero
    double tmp_element;
    for(i = 0; i< k*n; ++i) {
      if ( (*h0)[i] == 0. || numerh[i]  == 0.)
	h[i] = 0.;
      else {
	tmp_element = (*h0)[i] * (numerh[i] / (work2[i] + DIV_BY_ZERO_AVOIDANCE));
	h[i] = (tmp_element < ZERO_THRESHOLD) ? 0. : tmp_element;
      }
    }
    


    // calculating matrix w
    //----------------------------
    // calculating numerw = a*h'
    dgemm('N', 'T', m, k, n, 1.0, a, m, h, k, 0., numerw, m);
    // calculating first intermediate result work1 = h*h' (kxk-Matrix) => re-use of work1
    dgemm('N', 'T', k, k, n, 1.0, h, k, h, k, 0., work1, k);
    // calculating second intermediate result work2w = w0 * work1
    dgemm('N', 'N', m, k, k, 1.0, *w0, m, work1, k, 0., work2w, m);




    //calculating elementwise matrixmultiplication, Division and addition w = w0 .* (numerw ./ (work2w + eps))
    //set elements < zero_threshold to zero
    for(i = 0; i < m*k; ++i) {
      if ( (*w0)[i] == 0. || numerw[i] == 0.)
	w[i] = 0.;
      else {
	tmp_element = (*w0)[i] * (numerw[i] / (work2w[i] + DIV_BY_ZERO_AVOIDANCE));
	w[i] = (tmp_element < ZERO_THRESHOLD) ? 0. : tmp_element;
      }
    }  



    
    // calculating the norm of D = A-W*H
    dnorm = calculateNorm(a, w, h, d, m, n, k);

    
    // calculating change in w -> dw
    //----------------------------------
    double dw = calculateMaxchange(w, *w0, m, k, sqrteps);

    
    // calculating change in h -> dh
    //-----------------------------------
    double dh = calculateMaxchange(h, *h0, k, n, sqrteps);

    //Max-Change = max(dh, dw) = delta
    double delta = 0.0;
    delta = (dh > dw) ? dh : dw;


    // storing the matrix results of the current iteration in W0 respectively H0
    swap(w0, &w);
    swap(h0, &h);





    //Check for Convergence
    if (iter > 1) {
      if (delta < TolX) {
	*maxiter = iter;
	break;
      }
      else
	if (dnorm <= TolFun*dnorm0) {
	*maxiter = iter;
	break;
	}


    }

    
    // storing the norm result of the current iteration
    dnorm0 = dnorm;






#if DEBUG_LEVEL >= 1
  printf("iter: %.6d\t dnorm: %.16f\t delta: %.16f\n", iter, dnorm, delta);
#endif   
   
  } //end of loop from 1 to maxiter

#if DEBUG_LEVEL >= 2
	printf("Exiting nmf_mu\n");
#endif
#ifdef PROFILE_NMF_MU
	gettimeofday(&end, 0);
	outputTiming("", start, end);
#endif

  // freeing memory if used
  
  free(numerh);
  free(numerw);
  free(work1);
  free(work2);
  free(h);
  free(w);
  free(work2w);
  free(d);

  // returning calculated norm
  return dnorm;
}
//end of nmf_mu
//-------------
