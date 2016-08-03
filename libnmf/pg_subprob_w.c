/**
 * PG_SUBPROB_W	calculates the projected gradient step of factor matrix W used in nmf_pg and nmf_alspg
 * 
 * Purpose:
 *		This routine calculates the projected gradient step of factor matrix W used in nmf_pg and nmf_alspg
 *		It is split in two routines for W and H respectively to avoid unecessary transposition of input and output matrices
 *
 * Description:
 *		The general gradient descent approach is shown in following code:
 *
 *		W = rand(m,k);
 *		H = rand(k,n);
 *		for iter=1:maxiter
 *			H = H - epsH * gradH;
 *			W = W - epsW * gradW;
 *		end
 *
 *		This routine computes the necessary stepsize epsW. Using this stepsize the next iteration result
 *		of W is calculated.
 *		
 *		The algorithm is described in Lin (2007) (see references [2]).
 *
 *		This routine uses the BLAS "dgemm" routine for matrix-matrix multiplications
 *
 * Arguments:
 *
 * a		in	pointer to matrix to factorise
 *
 * w0		in/out	pointer to factor matrix w
 *			on exit w0 contains the new solution to w
 *
 * h		in	pointer to factor matrix h
 *
 * grad		in/out	pointer to allocated memory for the gradient
 *			on exit grad contains the gradient
 *
 * wn		in	pointer to allocated memory for possible next w
 *
 * aht		in	pointer to allocated memory for A*H'
 *
 * hht		in	pointer to allocated memory for H*H'
 *
 * d		in	pointer to allocated memory for D = A - W*H
 *
 * tempw	in	pointer to allocated memory for a temporary variable
 *
 * m		in, 	first dimension of a and w/w0
 *
 * n		in, 	second dimension of a and h/h0
 *
 * k		in, 	second dimension of w/w0 and first dimension of h/h0
 *
 * tol		in	tolerance value for stopping condition
 *
 * maxiter	in, 	maximal number of iterations to run
 *
 * function value	out, 	number of iterations actually run
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



int pg_subprob_w(double *a, double *w0, double *h, double *grad, double *wn, double *aht, double *hht, double * d, double * tempw, int m, int n, int k, double tol, int maxiter) 
{
  

  double * wp = (double*) malloc(sizeof(double)*m*k);
  int decr_factor_a;

  //aht = A * H'
  dgemm('N', 'T', m, k, n, 1.0, a, m, h, k, 0., aht, m);

  //hht = H * H'
  dgemm('N', 'T', k, k, n, 1.0, h, k, h, k, 0., hht, k);

  double factor_a = 1;
  double factor_b = 0.1;
  //wp == wn?
  int wpeqwn = 0;

  //Loop-Indices
  int i,iter, inner_iter;
  for (iter=1; iter<= maxiter; ++iter) {
    //grad = W * HHt - AHt
    dlacpy('A', m, k, aht, m, grad, m);
    dgemm('N', 'N', m, k, k, 1.0, w0, m, hht, k, -1.0, grad, m);
    

    double projgrad = 0;
    for(i = 0; i<m*k; ++i)
      if (grad[i] < 0 || w0[i] > 0)
	projgrad += grad[i]*grad[i];
    projgrad = sqrt(projgrad);

    //stopping criterion
    if (projgrad < tol)
      break;



    //Search step size
    //----------------
    for (inner_iter = 1; inner_iter <= 20; ++inner_iter) {

      //wn = max(W - factor_a * grad, 0);
      //hn = max(H-factor_a*grad, 0)
      
      // copy w0 to wn
      dlacpy('A', m, k, w0, m, wn, m);
      // compute wn = - factor_a * grad + wn
      daxpy(m*k, - factor_a, grad, 1, wn, 1);
      //set small positive and negative elements to zero for performance reasons and to keep the non-negativity constraint
      for (i = 0; i < m*k; ++i) {
	if (wn[i] < ZERO_THRESHOLD)
	  wn[i] = 0.;
      }
   

      //calculate difference between w of current and last iter

      //copy wn to d
      dlacpy('A', m, k, wn, m, d, m);



      // d = wn - w0
      daxpy(m*k, - 1.0, w0, 1, d, 1);
      
      // gradd = sum(sum(grad.*d));
      double gradd = 0;
      for (i = 0; i < m*k; ++i) {
	gradd += grad[i] * d[i];
      }
      
      
      //dqd = sum(sum((hht*d).*d))
      //tempw = d * hht = m x k
      dgemm('N', 'N', m, k, k, 1.0, d, m, hht, k, 0., tempw, m);
      //sum up elements of tempw .* d
      double dqd = 0.;
      for (i=0; i< m*k; ++i) {
	dqd += tempw[i]*d[i];
      }

      // test if sufficient decrease?!
      int suff_decr = (0.99 * gradd + 0.5*dqd) < 0;


      if (inner_iter == 1) {
	decr_factor_a = !suff_decr;
	//wp = w0;
	dlacpy('A', m, k, w0, m, wp, m);
	
      }

      if (decr_factor_a) {
	if (suff_decr) {
	  //w = wn;
	  dlacpy('A', m, k, wn, m, w0, m);
	  break;
	}
	else
	  factor_a = factor_a * factor_b;
      }
      else {
	//check wp == wn?
	wpeqwn = 1;
	for (i = 0; i < m*k; ++i) {
	  if (wp[i] != wn[i]) {
	    wpeqwn = 0;
	    break;
	  }
	}

	if (!suff_decr || wpeqwn) {
	  //w = wp;
	  dlacpy('A', m, k, wp, m, w0, m);
	  break;
	}
	else {
	  factor_a = factor_a/factor_b;
	  //wp = wn;
	  dlacpy('A', m, k, wn, m, wp, m);
	}
      }
    }
  }

  free(wp);
  return iter;
}
//end of pg_subprob
//-----------------
