/**
 * NMF_PG -	calculates the nmf using a direct projected gradient approach
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
 *		This routine implements an direct projected gradients gradients approach.
 *		Its algorithm is explained in Lin (2007) (see references [2])
 *
 *		The general gradient descent approach is shown in following code:
 *
 *		W = rand(m,k);
 *		H = rand(k,n);
 *		for iter=1:maxiter
 *			H = H - epsH * gradH;
 *			W = W - epsW * gradW;
 *		end
 *		
 *		The gradients are calculated as follows:
 *
 *		gradW = W * (H*H') - A*H';
 *		gradH = (W'*W)*H - W'*V;
 *
 *		This routine uses following BLAS/LAPACK routines:
 *		
 *		dgemm	matrix-matrix multiplications to calculate gradients
 *		dlange	used to compute the frobenius norm
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
 * tol		in, 	used in check for convergence, relative stopping criterion
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
#include "outputtiming.h"
#include "calculatenorm_singleprec.h"
#include "calculatemaxchange_singleprec.h"
#include "pg_subprob_h_singleprec.h"










float nmf_pg_singleprec(float ** a, float ** w0, float **h0, int m, int n, int k, int * maxiter, const float tol)
{
  

#ifdef PROFILE_NMF_PG
	struct timeval start, end;
	gettimeofday(&start, 0);
#endif
#if DEBUG_LEVEL >= 2
	printf("Entering nmf_pg\n");
#endif


#ifdef ERROR_CHECKING
  errno = 0;
#endif


  float * gradw = (float*) malloc(sizeof(float)*m*k);
  float * gradh = (float*) malloc(sizeof(float)*k*n);
  float * w = (float*) malloc(sizeof(float)*m*k);
  float * h = (float*) malloc(sizeof(float)*k*n);
  float * help1 = (float*) malloc(sizeof(float)*k*k);
  float * helph = (float*) malloc(sizeof(float)*k*n);
  float * helpw = (float*) malloc(sizeof(float)*k*m);
  float * hp = (float*) malloc(sizeof(float)*k*n);
  float * hn = (float*) malloc(sizeof(float)*k*n);
  float * wn = (float*) malloc(sizeof(float)*m*k);
  float * tempw = (float*) malloc(sizeof(float)*k*m);
  float * temph = (float*) malloc(sizeof(float)*k*n);

    //...for calculating the norm of A-W*H
  float* d = (float*) malloc(sizeof(float)*m*n);					//d = a - w*h
  float dnorm0 = 0;
  float dnorm = 0;
  

  //...for calculating maxChange
  float* diffh = (float*) malloc(sizeof(float)*k*n);
  //-------------------
  
  
  //Loop indices
  int iter, i;

#ifdef ERROR_CHECKING
  if (errno) {
    perror("Error allocating memory in nmf_pg");
    free(gradw);
    free(gradh);
    free(w);
    free(h);
    free(help1);
    free(helph);
    free(helpw);
    free(hp);
    free(hn);
    free(wn);
    free(tempw);
    free(temph);
    free(d);
    free(diffh);
    return -1;
  }
#endif
  
  float factor_a = 1;
  float initgrad = 0;
  float obj, newobj;

  //copying initialised starting matrices to w and h
  slacpy('A', m, k, *w0, m, w, m);

  slacpy('A', k, n, *h0, k, h, k);



#ifndef ZERO_THRESHOLD
  const float ZERO_THRESHOLD = eps;
#endif


  int iterw, iterh;

  // Begin of Iterations
  for (iter = 1; iter <= *maxiter; ++iter) {
  


  // calculate gradw = w*(h*h')-a*h'
  //--------------------------------

    //help1 = h*h'
    sgemm('N', 'T', k, k, n, 1.0, h, k, h, k, 0., help1, k);
    //gradw = w * help1
    sgemm('N', 'N', m, k, k, 1.0, w, m, help1, k, 0., gradw, m);
    //gradw = -a*h' + gradw
    sgemm('N', 'T', m, k, n, -1.0, *a, m, h, k, 1.0, gradw, m);


    // calculate gradh = (w'*w)*h-w'*a
    //--------------------------------

    //help1 = w'*w
    sgemm('T', 'N', k, k, m, 1.0, w, m, w, m, 0., help1, k);
    //gradh = help1 * h
    sgemm('N', 'N', k, n, k, 1.0, help1, k, h, k, 0., gradh, k);
    //gradh = -w'*a + gradh
    sgemm('T', 'N', k, n, m, -1.0, w, m, *a, m, 1.0, gradh, k);
    
    
    if (iter == 1) {
      
      initgrad = 0;
      initgrad += pow(slange('F', m, k, gradw, m, NULL), 2);
      initgrad += pow(slange('F', k, n, gradh, k, NULL), 2);
      initgrad = sqrt(initgrad);

      // H = pg_subprob_h
      swap_singleprec(h0, &h);
      iterh = pg_subprob_h_singleprec(*a, w, *h0, diffh, h, hp, help1, helph, temph, m,n,k, 0.001, 1000);

      
      dnorm = calculateNorm_singleprec(*a, w, h, d, m, n, k);
      obj = pow(dnorm, 2)*m*n/2;


      
#if DEBUG_LEVEL >= 1
  printf("iter: %.6d\t dnorm: %.16f\n", iter, dnorm);
#endif   

      continue;
    }
    
    
    float projnorm = 0; 

    // norm([gradW(gradW <0 | W > 0); gradH(gradH <0 | H > 0)]);
    for(i = 0; i<m*k; ++i)
      if (gradw[i] < 0 || (w)[i] > 0)
	projnorm += gradw[i]*gradw[i];
    for(i = 0; i<k*n; ++i)
      if (gradh[i] < 0 || (h)[i] > 0)
	projnorm += gradh[i]*gradh[i];
    projnorm = sqrt(projnorm);
    

    // stopping condition
    if (projnorm < (tol*initgrad)) {
      break;
    }


   
    //wn = max(W-factor_a*gradw, 0)
    
    // copy w to wn
    slacpy('A', m, k, w, m, wn, m);

    //calculate wn = - factor_a * gradw + wn
    saxpy(k*m, -factor_a, gradw, 1, wn, 1);
    
    // set negative matrix elements to zero
    for (i = 0; i < m*k; ++i) {
      if (wn[i] < 0.)
	wn[i] = 0.;
    }


    //hn = max(H-factor_a*gradh, 0)

    // copy h to hn
    slacpy('A', k, n, h, k, hn, k);

    //calculate hn = -factor_a * gradh + hn
    saxpy(k*n, -factor_a, gradh, 1, hn, 1);

    // set negative matrix elements to zero
    for (i = 0; i < k*n; ++i) {
      if (hn[i] < ZERO_THRESHOLD)
	hn[i] = 0.;
    }


    //adjust factor_a
    //----------------

    newobj = pow(calculateNorm_singleprec(*a, wn, hn, d, m, n, k), 2)*m*n/2;
    
    
    float compval = 0;
    for (i = 0; i < k*m; ++i)
      compval += gradw[i] * (wn[i] - w[i]);
    for (i = 0; i < k*n; ++i)
      compval += gradh[i] * (hn[i] - h[i]);

    if ((newobj - obj) > 0.01 * compval ){	//decrease factor_a
      while(1) {
	factor_a /= 10;
	
	//wn = max(W-factor_a*gradw, 0)
	
	// copy w to wn
	slacpy('A', m, k, w, m, wn, m);

	//calculate wn = - factor_a * gradw + wn
	saxpy(k*m, -factor_a, gradw, 1, wn, 1);
	
	// set negative matrix elements to zero
	for (i = 0; i < m*k; ++i) {
	  if (wn[i] < ZERO_THRESHOLD)
	    wn[i] = 0.;
	}


	//hn = max(H-factor_a*gradh, 0)

	// copy h to hn
	slacpy('A', k, n, h, k, hn, k);

	//calculate hn = -factor_a * gradh + hn
	saxpy(k*n, -factor_a, gradh, 1, hn, 1);

	// set negative matrix elements to zero
	for (i = 0; i < k*n; ++i) {
	  if (hn[i] < ZERO_THRESHOLD)
	    hn[i] = 0.;
	}
	  
	newobj = pow(calculateNorm_singleprec(*a, wn, hn, d, m, n, k), 2)*m*n/2;

	compval = 0;
	for (i = 0; i < k*m; ++i)
	  compval += gradw[i] * (wn[i] - w[i]);
	for (i = 0; i < k*n; ++i)
	  compval += gradh[i] * (hn[i] - h[i]);
	
	if ((newobj - obj) <= 0.01*compval ){
	  swap_singleprec(&w, &wn);
	  swap_singleprec(&h, &hn);
	  obj = newobj;
	  break;
	}
      }
    }
    else {					//increase factor_a
      while(1) {
	//storing values in case they are the next iteration result
	slacpy('A', m, k, wn, m, helpw, m);
	
	slacpy('A', k, n, hn, k, helph, k);

	float objp = newobj;

	factor_a *= 10;

	
	//wn = max(W-factor_a*gradw, 0)
	
	// copy w to wn
	slacpy('A', m, k, w, m, wn, m);

	//calculate wn = - factor_a * gradw + wn
	saxpy(k*m, -factor_a, gradw, 1, wn, 1);
	
	// set negative matrix elements to zero
	for (i = 0; i < m*k; ++i) {
	  if (wn[i] < ZERO_THRESHOLD)
	    wn[i] = 0.;
	}


	//hn = max(H-factor_a*gradh, 0)

	// copy h to hn
	slacpy('A', k, n, h, k, hn, k);

	//calculate hn = -factor_a * gradh + hn
	saxpy(k*n, -factor_a, gradh, 1, hn, 1);

	// set negative matrix elements to zero
	for (i = 0; i < k*n; ++i) {
	  if (hn[i] < ZERO_THRESHOLD)
	    hn[i] = 0.;
	}

	newobj = pow(calculateNorm_singleprec(*a, wn, hn, d, m, n, k), 2)*m*n/2;

	//Value to compare to
	compval = 0;
	for (i = 0; i < k*m; ++i)
	  compval += gradw[i] * (wn[i] - w[i]);
	for (i = 0; i < k*n; ++i)
	  compval += gradh[i] * (hn[i] - h[i]);

	//Alternative condition [Wn Hn'] == [Wp Hp']
	int equal = 1;
	for (i = 0; i < m*k; ++i) {
	  if (wn[i] != helpw[i]) {
	    equal = 0;
	    break;
	  }
	}
	if (equal)
	  for(i = 0; i < k*n; ++i)
	    if (hn[i] != helph[i]) {
	      equal = 0;
	      break;
	    }

	if ( ((newobj - obj) > 0.01*compval) || equal) {
	  swap_singleprec(&w, &helpw);
	  swap_singleprec(&h, &helph);

	  obj = objp;
	  factor_a /= 10;
	  
	  break;
	}
      }




    }



  } //end of for

  
    // calculating the norm of D = A-W*H
    dnorm = calculateNorm_singleprec(*a, w, h, d, m, n, k);


    // storing the norm results of the current iteration
    dnorm0 = dnorm;
    

#if DEBUG_LEVEL >= 1
  printf("iter: %.6d\t dnorm: %.16f\t proj-grad norm: %.16f\n", iter, dnorm, projnorm);
#endif   
   


  // storing the matrix results of the current iteration
  swap_singleprec(w0, &w);
  swap_singleprec(h0, &h);


#if DEBUG_LEVEL >= 2
	printf("Exiting nmf_pg\n");
#endif
#ifdef PROFILE_NMF_PG
	gettimeofday(&end, 0);
	outputTiming("", start, end);
#endif

  // freeing memory if used
  free(gradw);
  free(gradh);
  free(w);
  free(h);
  free(help1);
  free(helph);
  free(helpw);
  free(hp);
  free(hn);
  free(wn);
  free(tempw);
  free(temph);
  free(d);
  free(diffh);


  *maxiter = iter;
  // returning calculated norm
  return dnorm;
}
//end of nmf_pg
//-------------
