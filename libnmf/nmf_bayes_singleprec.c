/**
 * NMF_BAYES -	calculates the nmf using a bayesion priors approach
 *
 *
 
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
 *		Algorithm by Mikkel Schmidt, see Documentation for References
 *			
 * Arguments:
 *
 * a		in, 	pointer to matrix to factorise
 *
 * w0		in, 	pointer to initialised factor-matrix w0
 *		out, 	pointer to final factor-matrix w
 * h0		in, 	pointer to initialised factor-matrix h0
 *		out, 	pointer to final factor-matrix h
 * m		in, 	first dimension of a and w/w0
 *
 * n		in, 	second dimension of a and h/h0
 *
 * k		in, 	second dimension of w/w0 and first dimension of h/h0
 *
 * maxiter	in, 	maximal number of iterations to run
 *		out, 	number of iterations actually run
 *
 * TolX		in, 	used in check for convergence, tolerance for maxchange
 *
 * TolFun	in, 	used in check for convergence, tolerance for dnorm
 *
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
#include "randnumber_singleprec.h"
#include "storematrix_singleprec.h"









float nmf_bayes_singleprec(float ** a, float ** w0, float ** h0, int m, int n, \
		      int k, int * maxiter, const float TolX, const float TolFun) 
{
  

#ifdef PROFILE_NMF_NEALS
	struct timeval start, end;
	gettimeofday(&start, 0);
#endif
#if DEBUG_LEVEL >= 2
	printf("Entering nmf_neals\n");
#endif


#ifdef ERROR_CHECKING
  errno = 0;
#endif


  float * help1 = (float*) malloc(sizeof(float)*k*k);
  float * help2 = (float*) malloc(sizeof(float)*k*n);
  float * help3 = (float*) malloc(sizeof(float)*k*m);
  float * help4 = (float*) malloc(sizeof(float)*m*(k-1)); //W(m,k-1)
  float * help5 = (float*) malloc(sizeof(float)*(k-1));   //(B*B')(k-1,1)
  float * help6 = (float*) malloc(sizeof(float)*(m));     //col vector of W/alpha	
  float * help7 = (float*) malloc(sizeof(float)*(k-1)*n); //H(k-1,n)
  float * help8 = (float*) malloc(sizeof(float)*n);	     //row vector 1xn
  int * zerocolrow = (int*) malloc(sizeof(float) * k);	     //indicates if a whole row/column is zero (1) or not (0)
  //-----------------------------------------



  // definition of necessary dynamic data structures
  //...for calculating matrix h
  float* h = (float*) malloc(sizeof(float)*k*n);
  float* beta = (float*) malloc(sizeof(float)*k*n);
  slaset('A', k, n, 0.0, 0.0, beta, k);
  
  //...for calculating matrix w
   float* w = (float*) malloc(sizeof(float)*m*k);
   float* alpha = (float*) malloc(sizeof(float)*m*k);
   slaset('A', m, k, 0.0, 0.0, alpha, m);
  //----------------

  // -- Werte einlesen
  int kfactor = 0;
  int maxiiter = 20;
  float theta = 0.;
  float sigma = 1.;
  

  //...for calculating the norm of A-W*H
  float* d = (float*) malloc(sizeof(float)*m*n);					//d = a - w*h
  
  float dnorm = 0;
  

  //-------------------

#ifdef ERROR_CHECKING
  if (errno) {
    perror("Error allocating memory in nmf_neals");
    free(help1);
    free(help2);
    free(help3);
    free(help4);
    free(help5);
    free(help6);
    free(help7);
    free(help8);
    free(zerocolrow);
    free(alpha);
    free(beta);
    free(h);
    free(w);
    free(d);
    return -1;
  }
#endif


   
  //Loop-Indices
  int iter, iiter, i;
  
  // x = sum(A(:).^2)/2 = frobnorm(A)^2 / 2
  float x = 0.;
  x = slange('F', m, n, *a, m, NULL);
  x *= x;
  x /= 2.;
  

   // factorisation step in a loop from 1 to maxiter
  for (iter = 1; iter <= *maxiter; ++iter) {



    // calculating matrix w update
    //----------------------------
    //help1 = h*h'
    sgemm('N', 'T', k, k, n, 1.0, *h0, k, *h0, k, 0., help1, k);
    //help3 = a*h'
    sgemm('N', 'T', m, k, n, 1.0, *a, m, *h0, k, 0., help3, m);

    
    // inner iterations
    for (iiter = 1; iiter <= maxiiter; ++iiter) {
        //for every n from 1 to k update the n-th column of W
	int nk;
	for (nk = 1; nk <= k; ++nk) {
	  //copy every coklumn except for the n-th of W to help4
	  int tmpi, col = 0;
	  for (tmpi = 1; tmpi <=k; ++tmpi) {
	    if (tmpi != nk) {
	      scopy(m, (*w0) + (tmpi-1) * m, 1, help4 + col * m, 1);
	      ++col;
	    }
	  }
	  //copy every row except for the n-th of column n of H*H' to help5
	  int row = 0;
	  for (tmpi = 1; tmpi <=k; ++tmpi) {
	    if (tmpi != nk) {
	      help5[row] = help1[(tmpi-1) + (nk-1) * k];
	      ++row;
	    }
	  }

	  // calculate help4*help5
	  //help6 = -help4*help5
	  sgemm('N', 'N', m, 1, k-1, -1.0, help4, m, help5, k-1, 0.0, help6, m);
	  // scale alpha(:,n)
	  //help6 = help6 - sigma * alpha(:,n) / (h*h')(n,n)
	  saxpy(m, -sigma, alpha + (nk-1) * m, 1, help6, 1);
	  // calculate help3 - help4*help5 - scaled alpha
	  
	  
	  int tmpj, wholecolzero = 1;
	  float newval;
	  for(tmpj = 0; tmpj < m; ++tmpj) {
	    newval = (help3[tmpj + (nk-1) * m] + help6[tmpj])/help1[(nk-1) + (nk-1) * k];
	    (*w0)[tmpj + (nk-1) * m] = newval > ZERO_THRESHOLD ? newval : 0.;
	    if (newval > ZERO_THRESHOLD)
	      wholecolzero = 0;
	  }
	  zerocolrow[nk-1] = wholecolzero;
	}
    }
 
    //set all columns of W which consist only of zeros to mean(alpha,2).*(-log(rand(m, {count of "zero columns"})))
    
    //create replacement vector help6 = mean(alpha,2)
    int tmpi, tmpj;
    float mean;
    for(tmpi = 0; tmpi < m; ++tmpi) {
      mean = 0.;
      for(tmpj = 0; tmpj < k; ++tmpj) {
	mean += alpha[tmpi + tmpj * m];
      }
      help6[tmpi] = mean / k;
    }
    // run over every column of W and...
    // check zerocolrow[currentColumn] to see if it is to be updated
    for(tmpi = 0; tmpi < k; ++tmpi) {
	if (zerocolrow[tmpi]) {
	    //set column to help6 .* (-log(rand(m,1)))
	    for(tmpj = 0; tmpj < k; ++tmpj) {
		(*w0)[tmpj + tmpi * m] = help6[tmpj] * (-log(randnumber_singleprec(0,1)));
	    }
	}
    }
    
 
    
    //sigma update
    // calculate help3 = W*(H*H')- 2 * (A*H') size m x k
    sgemm('N', 'N', m, k, k, 1.0, *w0, m, help1, k, -2.0, help3, m);
    //calculate sum( sum( W .* help3 ) )
    float sum = 0;
    for(tmpi = 0; tmpi < m*k; ++tmpi) {
	sum += help3[tmpi] * (*w0)[tmpi]; 
    }
    
    sigma = (theta + x + sum/2) / (m*n/2 + kfactor + 1);
   
    

    
    // calculating matrix h
    //----------------
    //help1 = w0'*w0
    sgemm('T', 'N', k, k, m, 1.0, *w0, m, *w0, m, 0., help1, k);
    //help2 = w0'*a
    sgemm('T', 'N', k, n, m, 1.0, *w0, m, *a, m, 0., help2, k);




    
    // inner iterations
    for (iiter = 1; iiter <= maxiiter; ++iiter) {
        //for every n from 1 to k update the n-th column of W
	int nk;
	for (nk = 1; nk <= k; ++nk) {
	  //copy every row except for the n-th of H to help7
	  int tmpi, row = 0;
	  for (tmpi = 1; tmpi <=k; ++tmpi) {
	    if (tmpi != nk) {
	      scopy(n, (*h0) + (tmpi-1), k, help7 + row, (k-1));
	      ++row;
	    }
	  }
	  //copy every col except for the n-th of row n of W'*W (help1) to help5
	  int col = 0;
	  for (tmpi = 1; tmpi <=k; ++tmpi) {
	    if (tmpi != nk) {
	      help5[col] = help1[(nk - 1) + (tmpi -1) * k];
	      ++col;
	    }
	  }
	  // calculate help5*help7
	  //help8 = -help5*help7 (size 1xn)
	  sgemm('N', 'N', 1, n, k-1, -1.0, help5, 1, help7, k-1, -0.0, help8, 1);
	  // scale beta(n,:)
	  //help8 = help8 - sigma * beta(n,:) / (w'*w)(n,n)
	  saxpy(n, -sigma, beta + nk-1, k, help8, 1);
	  // calculate help2 - help5*help7 - scaled beta
	  int tmpj, wholerowzero = 1;
	  float newval;
	  for(tmpj = 0; tmpj < n; ++tmpj) {
	    newval = (help2[(nk-1) + tmpj * k] + help8[tmpj])/help1[(nk-1) + (nk-1) * k];
	    (*h0)[(nk-1) + tmpj * k] = newval > ZERO_THRESHOLD ? newval : 0.;
	    if (newval > ZERO_THRESHOLD)
	      wholerowzero = 0;
	  }
	  zerocolrow[nk-1] = wholerowzero;
	}
    }
    
    //set all rows of H which consist only of zeros to mean(beta,1).*(-log(rand({count of "zero rows"},n )))
    
    //create replacement vector help8 = mean(beta,1)

    for(tmpi = 0; tmpi < n; ++tmpi) {
      mean = 0.;
      for(tmpj = 0; tmpj < k; ++tmpj) {
	mean += beta[tmpj + tmpi * k];
      }
      help8[tmpi] = mean / k;
    }
    // run over every row of H and...
    // check zerocolrow[currentRow] to see if it is to be updated
    for(tmpi = 0; tmpi < k; ++tmpi) {
	if (zerocolrow[tmpi]) {
	    //set row to help8 .* (-log(rand(m,1)))
	    for(tmpj = 0; tmpj < n; ++tmpj) {
		(*h0)[tmpi + tmpj * k] = help8[tmpj] * (-log(randnumber_singleprec(0,1)));
	    }
	}
    }      





   
  } //end of loop from 1 to maxiter

    
    // calculating the norm of D = A-W*H
    dnorm = calculateNorm_singleprec(*a, *w0, *h0, d, m, n, k);
    


    

#if DEBUG_LEVEL >= 1
  printf("iter: %.6d\t dnorm: %.16f\t delta: %.16f\n", iter, dnorm, delta);
#endif   

#if DEBUG_LEVEL >= 2
	printf("Exiting nmf_bayes\n");
#endif
#ifdef PROFILE_NMF_NEALS
	gettimeofday(&end, 0);
	outputTiming("", start, end);
#endif

  // freeing memory if used
    free(help1);
    free(help2);
    free(help3);
    free(help4);
    free(help5);
    free(help6);
    free(help7);
    free(help8);
    free(h);
    free(w);
    free(d);
    free(zerocolrow);
    free(alpha);
    free(beta);

    
  // returning calculated norm
  return dnorm;
}
//end of nmf_neals
//-------------
