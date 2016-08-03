/**
 * PG_SUBPROB_H calculates the projected gradient step of factor matrix H used in nmf_pg and nmf_alspg
 * 
 * Purpose:
 *		This routine calculates the projected gradient step of factor matrix H used in nmf_pg and nmf_alspg
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
 *		This routine computes the necessary stepsize epsH. Using this stepsize the next iteration result
 *		of H is calculated.
 *		
 *		The algorithm is described in Lin (2007) (see references [2]).
 *
 *		This routine uses the BLAS "dgemm" routine for matrix-matrix multiplications
 *
 * Arguments:
 *
 * a		in	pointer to matrix to factorise
 *
 * w		in	pointer to factor matrix w
 *
 * h0		in/out	pointer to factor matrix h
 *			on exit h0 contains the new solution to h
 *
 * grad		in/out	pointer to allocated memory for the gradient
 *			on exit grad contains the gradient
 *
 * hn		in	pointer to allocated memory for possible next h
 *
 * wta		in	pointer to allocated memory for W'*A
 *
 * wtw		in	pointer to allocated memory for W'*W
 *
 * d		in	pointer to allocated memory for D = A - W*H
 *
 * temph	in	pointer to allocated memory for a temporary variable
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



int pg_subprob_h(double *a, double *w, double *h0, double *grad, double *hn, double *wta, double *wtw, double * d, double * temph, int m, int n, int k, double tol, int maxiter) 
{
  

  double * hp = (double*) malloc(sizeof(double)*k*n);
  int decr_factor_a;

  //wta = W' * A;
  dgemm('T', 'N', k, n, m,  1.0, w, m, a, m, 0., wta, k);

  //wtw = w'*w;
  dgemm('T', 'N', k, k, m, 1.0, w, m, w, m, 0., wtw, k);

  double factor_a = 1;
  double factor_b = 0.1;
  //hp == hn?
  int hpeqhn = 0;

  //Loop-Indices
  int i,iter, inner_iter;
  for (iter=1; iter<= maxiter; ++iter) {
    //grad = wtw*h-wta;
    dlacpy('A', k, n, wta, k, grad, k);
    
    dgemm('N', 'N', k, n, k, 1.0, wtw, k, h0, k, - 1.0, grad, k);

    double projgrad = 0;
    for(i = 0; i<k*n; ++i)
      if (grad[i] < 0 || h0[i] > 0)
	projgrad += grad[i]*grad[i];
    projgrad = sqrt(projgrad);

    //stopping criterion
    if (projgrad < tol)
      break;

    //Search step size
    //----------------
    for (inner_iter = 1; inner_iter <= 20; ++inner_iter) {

      
      //hn = max(H-factor_a*grad, 0)
      
      // copy h0 to hn
      dlacpy('A', k, n, h0, k, hn, k);
      // computy hn = - factor_a * grad + hn
      daxpy(k*n, - factor_a, grad, 1, hn, 1);
      //set small positive and negative elements to zero for performance reasons and do keep the non-negativity constraintd
      for (i = 0; i < k*n; ++i) {
	if (hn[i] < ZERO_THRESHOLD)
	  hn[i] = 0.;
      }
   

      

      //copy hn to d
      dlacpy('A', k, n, hn, k, d, k);

      // d = hn - h0
      daxpy(k*n, - 1.0, h0, 1, d, 1);
      
      // gradd = sum(sum(grad.*d));
      double gradd = 0;
      for (i = 0; i < k*n; ++i) {
	gradd += grad[i] * d[i];
      }
      
      
      //dqd = sum(sum((wtw*d).*d))
      dgemm('N', 'N', k, n, k, 1.0, wtw, k, d, k, 0., temph, k);
      double dqd = 0;
      for (i=0; i< k*n; ++i) {
	dqd += temph[i]*d[i];
      }

      // test if sufficient decrease?!
      int suff_decr = (0.99 * gradd + 0.5*dqd) < 0;


      if (inner_iter == 1) {
	decr_factor_a = !suff_decr;
	//hp = h;
	dlacpy('A', k, n, h0, k, hp, k);
	
      }

      if (decr_factor_a) {
	if (suff_decr) {
	  //h = hn;
	  dlacpy('A', k, n, hn, k, h0, k);
	  break;
	}
	else
	  factor_a = factor_a * factor_b;
      }
      else {
	//check hp == hn?
	hpeqhn = 1;
	for (i = 0; i < k*n; ++i) {
	  if (hp[i] != hn[i]) {
	    hpeqhn = 0;
	    break;
	  }
	}

	if (!suff_decr || hpeqhn) {
	  //h = hp;
	  dlacpy('A', k, n, hp, k, h0, k);
	  break;
	}
	else {
	  factor_a = factor_a/factor_b;
	  //hp = hn;
	  dlacpy('A', k, n, hn, k, hp, k);
	}
      }
    }
  }

  free(hp);
  return iter;
}
//end of pg_subprob
//-----------------
