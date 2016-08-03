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
#include "calculatenorm.h"
#include "calculatemaxchange.h"
#include "pg_subprob_h.h"










double nmf_pg(double ** a, double ** w0, double **h0, int m, int n, int k, int * maxiter, const double tol)
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


  double * gradw = (double*) malloc(sizeof(double)*m*k);
  double * gradh = (double*) malloc(sizeof(double)*k*n);
  double * w = (double*) malloc(sizeof(double)*m*k);
  double * h = (double*) malloc(sizeof(double)*k*n);
  double * help1 = (double*) malloc(sizeof(double)*k*k);
  double * helph = (double*) malloc(sizeof(double)*k*n);
  double * helpw = (double*) malloc(sizeof(double)*k*m);
  double * hp = (double*) malloc(sizeof(double)*k*n);
  double * hn = (double*) malloc(sizeof(double)*k*n);
  double * wn = (double*) malloc(sizeof(double)*m*k);
  double * tempw = (double*) malloc(sizeof(double)*k*m);
  double * temph = (double*) malloc(sizeof(double)*k*n);

    //...for calculating the norm of A-W*H
  double* d = (double*) malloc(sizeof(double)*m*n);					//d = a - w*h
  double dnorm0 = 0;
  double dnorm = 0;
  

  //...for calculating maxChange
  double* diffh = (double*) malloc(sizeof(double)*k*n);
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
  
  double factor_a = 1;
  double initgrad = 0;
  double obj, newobj;

  //copying initialised starting matrices to w and h
  dlacpy('A', m, k, *w0, m, w, m);

  dlacpy('A', k, n, *h0, k, h, k);



#ifndef ZERO_THRESHOLD
  const double ZERO_THRESHOLD = eps;
#endif


  int iterw, iterh;

  // Begin of Iterations
  for (iter = 1; iter <= *maxiter; ++iter) {
  


  // calculate gradw = w*(h*h')-a*h'
  //--------------------------------

    //help1 = h*h'
    dgemm('N', 'T', k, k, n, 1.0, h, k, h, k, 0., help1, k);
    //gradw = w * help1
    dgemm('N', 'N', m, k, k, 1.0, w, m, help1, k, 0., gradw, m);
    //gradw = -a*h' + gradw
    dgemm('N', 'T', m, k, n, -1.0, *a, m, h, k, 1.0, gradw, m);


    // calculate gradh = (w'*w)*h-w'*a
    //--------------------------------

    //help1 = w'*w
    dgemm('T', 'N', k, k, m, 1.0, w, m, w, m, 0., help1, k);
    //gradh = help1 * h
    dgemm('N', 'N', k, n, k, 1.0, help1, k, h, k, 0., gradh, k);
    //gradh = -w'*a + gradh
    dgemm('T', 'N', k, n, m, -1.0, w, m, *a, m, 1.0, gradh, k);
    
    
    if (iter == 1) {
      
      initgrad = 0;
      initgrad += pow(dlange('F', m, k, gradw, m, NULL), 2);
      initgrad += pow(dlange('F', k, n, gradh, k, NULL), 2);
      initgrad = sqrt(initgrad);

      // H = pg_subprob_h
      swap(h0, &h);
      iterh = pg_subprob_h(*a, w, *h0, diffh, h, hp, help1, helph, temph, m,n,k, 0.001, 1000);

      
      dnorm = calculateNorm(*a, w, h, d, m, n, k);
      obj = pow(dnorm, 2)*m*n/2;


      
#if DEBUG_LEVEL >= 1
  printf("iter: %.6d\t dnorm: %.16f\n", iter, dnorm);
#endif   

      continue;
    }
    
    
    double projnorm = 0; 

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
    dlacpy('A', m, k, w, m, wn, m);

    //calculate wn = - factor_a * gradw + wn
    daxpy(k*m, -factor_a, gradw, 1, wn, 1);
    
    // set negative matrix elements to zero
    for (i = 0; i < m*k; ++i) {
      if (wn[i] < 0.)
	wn[i] = 0.;
    }


    //hn = max(H-factor_a*gradh, 0)

    // copy h to hn
    dlacpy('A', k, n, h, k, hn, k);

    //calculate hn = -factor_a * gradh + hn
    daxpy(k*n, -factor_a, gradh, 1, hn, 1);

    // set negative matrix elements to zero
    for (i = 0; i < k*n; ++i) {
      if (hn[i] < ZERO_THRESHOLD)
	hn[i] = 0.;
    }


    //adjust factor_a
    //----------------

    newobj = pow(calculateNorm(*a, wn, hn, d, m, n, k), 2)*m*n/2;
    
    
    double compval = 0;
    for (i = 0; i < k*m; ++i)
      compval += gradw[i] * (wn[i] - w[i]);
    for (i = 0; i < k*n; ++i)
      compval += gradh[i] * (hn[i] - h[i]);

    if ((newobj - obj) > 0.01 * compval ){	//decrease factor_a
      while(1) {
	factor_a /= 10;
	
	//wn = max(W-factor_a*gradw, 0)
	
	// copy w to wn
	dlacpy('A', m, k, w, m, wn, m);

	//calculate wn = - factor_a * gradw + wn
	daxpy(k*m, -factor_a, gradw, 1, wn, 1);
	
	// set negative matrix elements to zero
	for (i = 0; i < m*k; ++i) {
	  if (wn[i] < ZERO_THRESHOLD)
	    wn[i] = 0.;
	}


	//hn = max(H-factor_a*gradh, 0)

	// copy h to hn
	dlacpy('A', k, n, h, k, hn, k);

	//calculate hn = -factor_a * gradh + hn
	daxpy(k*n, -factor_a, gradh, 1, hn, 1);

	// set negative matrix elements to zero
	for (i = 0; i < k*n; ++i) {
	  if (hn[i] < ZERO_THRESHOLD)
	    hn[i] = 0.;
	}
	  
	newobj = pow(calculateNorm(*a, wn, hn, d, m, n, k), 2)*m*n/2;

	compval = 0;
	for (i = 0; i < k*m; ++i)
	  compval += gradw[i] * (wn[i] - w[i]);
	for (i = 0; i < k*n; ++i)
	  compval += gradh[i] * (hn[i] - h[i]);
	
	if ((newobj - obj) <= 0.01*compval ){
	  swap(&w, &wn);
	  swap(&h, &hn);
	  obj = newobj;
	  break;
	}
      }
    }
    else {					//increase factor_a
      while(1) {
	//storing values in case they are the next iteration result
	dlacpy('A', m, k, wn, m, helpw, m);
	
	dlacpy('A', k, n, hn, k, helph, k);

	double objp = newobj;

	factor_a *= 10;

	
	//wn = max(W-factor_a*gradw, 0)
	
	// copy w to wn
	dlacpy('A', m, k, w, m, wn, m);

	//calculate wn = - factor_a * gradw + wn
	daxpy(k*m, -factor_a, gradw, 1, wn, 1);
	
	// set negative matrix elements to zero
	for (i = 0; i < m*k; ++i) {
	  if (wn[i] < ZERO_THRESHOLD)
	    wn[i] = 0.;
	}


	//hn = max(H-factor_a*gradh, 0)

	// copy h to hn
	dlacpy('A', k, n, h, k, hn, k);

	//calculate hn = -factor_a * gradh + hn
	daxpy(k*n, -factor_a, gradh, 1, hn, 1);

	// set negative matrix elements to zero
	for (i = 0; i < k*n; ++i) {
	  if (hn[i] < ZERO_THRESHOLD)
	    hn[i] = 0.;
	}

	newobj = pow(calculateNorm(*a, wn, hn, d, m, n, k), 2)*m*n/2;

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
	  swap(&w, &helpw);
	  swap(&h, &helph);

	  obj = objp;
	  factor_a /= 10;
	  
	  break;
	}
      }




    }



  } //end of for

  
    // calculating the norm of D = A-W*H
    dnorm = calculateNorm(*a, w, h, d, m, n, k);


    // storing the norm results of the current iteration
    dnorm0 = dnorm;
    

#if DEBUG_LEVEL >= 1
  printf("iter: %.6d\t dnorm: %.16f\t proj-grad norm: %.16f\n", iter, dnorm, projnorm);
#endif   
   


  // storing the matrix results of the current iteration
  swap(w0, &w);
  swap(h0, &h);


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
