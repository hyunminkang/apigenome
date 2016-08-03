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
 *		This routine uses the BLAS "sgemm" routine for matrix-matrix-multiplications
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
#include "calculatenorm_singleprec.h"
#include "calculatemaxchange_singleprec.h"




float nmf_mu_singleprec(float * a, float ** w0, float ** h0, int m, int n, \
		      int k, int * maxiter, const float TolX, const float TolFun) 
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
  float* numerh = (float*) malloc(sizeof(float) *k*n);
  float* work1 = (float*) malloc(sizeof(float)*k*k);					// used for calculation of h & w
  float* work2 = (float*) malloc(sizeof(float)*k*n);
  float* h = (float*) malloc(sizeof(float)*k*n);
  //----------------

  //...for calculating matrix w
  float* numerw = (float*) malloc(sizeof(float)*m*k);
  float* work2w = (float*) malloc(sizeof(float)*m*k);
  float* w = (float*) malloc(sizeof(float)*m*k);
  //-----------------

  //...for calculating the norm of A-W*H
  float* d = (float*) malloc(sizeof(float)*m*n);					//d = a - w*h
  float dnorm0 = 0;
  float dnorm = 0;
  const float eps = slamch('E');					//machine precision epsilon
  const float sqrteps = sqrt(eps);					//squareroot of epsilon
  


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
  const float ZERO_THRESHOLD = eps;
#endif

  
  //Loop-Indices
  int iter, i;

  // factorisation step in a loop from 1 to maxiter
  for (iter = 1; iter <= *maxiter; ++iter) {




    // calculating matrix h
    //----------------
    // calculating numerh = w0'*a
    sgemm('T', 'N', k, n, m, 1.0, *w0, m, a, m, 0., numerh, k);  
    // calculating first intermediate result work1 = w0'*w0
    sgemm('T', 'N', k, k, m, 1.0, *w0, m, *w0, m, 0., work1, k);
    // calculating second intermediate result work2 = work1 * h0
    sgemm('N', 'N', k, n, k, 1.0, work1, k, *h0, k, 0., work2, k);



    //calculating elementwise matrixmultiplication, Division and addition h = h0 .* (numerh ./(work2 + eps))
    //set elements < zero_threshold to zero
    float tmp_element;
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
    sgemm('N', 'T', m, k, n, 1.0, a, m, h, k, 0., numerw, m);
    // calculating first intermediate result work1 = h*h' (kxk-Matrix) => re-use of work1
    sgemm('N', 'T', k, k, n, 1.0, h, k, h, k, 0., work1, k);
    // calculating second intermediate result work2w = w0 * work1
    sgemm('N', 'N', m, k, k, 1.0, *w0, m, work1, k, 0., work2w, m);




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
    dnorm = calculateNorm_singleprec(a, w, h, d, m, n, k);

    
    // calculating change in w -> dw
    //----------------------------------
    float dw = calculateMaxchange_singleprec(w, *w0, m, k, sqrteps);

    
    // calculating change in h -> dh
    //-----------------------------------
    float dh = calculateMaxchange_singleprec(h, *h0, k, n, sqrteps);

    //Max-Change = max(dh, dw) = delta
    float delta = (dh > dw) ? dh : dw;


    // storing the matrix results of the current iteration in W0 respectively H0
    swap_singleprec(w0, &w);
    swap_singleprec(h0, &h);



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
