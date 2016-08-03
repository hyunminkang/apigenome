/**
 * NMF_ALSPG -	calculates the nmf using an alternating least square approach applying a project gradient method
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
 *		This routine implements an ALS approach using projected gradients to solve the subproblems of this
 *		approach. Its algorithm is explained in Lin (2007) (see references [2])
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
#include "pg_subprob_w.h"









double nmf_alspg(double ** a, double ** w0, double **h0, int m, int n, int k, int * maxiter, const double tol)
{
  

#ifdef PROFILE_NMF_ALSPG
        struct timeval start, end;
        gettimeofday(&start, 0);
#endif
#if DEBUG_LEVEL >= 2
	printf("Entering nmf_alspg\n");
#endif


#ifdef ERROR_CHECKING
  errno = 0;
#endif

  


  double * gradw = (double*) malloc(sizeof(double)*m*k);
  double * gradh = (double*) malloc(sizeof(double)*k*n);
  double * h = (double*) malloc(sizeof(double) * k*n);
  double * help1 = (double*) malloc(sizeof(double)*k*k);
  double * helph = (double*) malloc(sizeof(double)*k*n);
  double * helpw = (double*) malloc(sizeof(double)*k*m);
  double * hp = (double*) malloc(sizeof(double)*k*n);
  double * hn = (double*) malloc(sizeof(double)*k*n);
  double * hpw = (double*) malloc(sizeof(double)*k*m);
  double * hnw = (double*) malloc(sizeof(double)*k*m);
  double * tempw = (double*) malloc(sizeof(double)*k*m);
  double * temph = (double*) malloc(sizeof(double)*k*n);

    //...for calculating the norm of A-W*H
  double* d = (double*) malloc(sizeof(double)*m*n);					//d = a - w*h
  double dnorm0 = 0;
  double dnorm = 0;
  

  //...for calculating maxChange
  double* diffw = (double*) malloc(sizeof(double)*m*k);
  double* diffh = (double*) malloc(sizeof(double)*k*n);
  //-------------------

  
  
  //Loop indices
  int iter, i;

#ifdef ERROR_CHECKING
  if (errno) {
    perror("Error allocating memory in nmf_alspg");
    free(gradw);
    free(gradh);
    free(help1);
    free(helph);
    free(helpw);
    free(h);
    free(hp);
    free(hn);
    free(hpw);
    free(hnw);
    free(tempw);
    free(temph);
    free(d);
    free(diffw);
    free(diffh);
    return -1;
  }
#endif
  
  // calculate gradw = w*(h*h')-a*h'
  //--------------------------------

  //help1 = h*h'
  dgemm('N', 'T', k, k, n, 1.0, *h0, k, *h0, k, 0., help1, k);
  //gradw = w * help1
  dgemm('N', 'N', m, k, k, 1.0, *w0, m, help1, k, 0., gradw, m);
  //gradw = -a*h' + gradw
  dgemm('N', 'T', m, k, n, -1.0, *a, m, *h0, k, 1.0, gradw, m);

  // calculate gradh = (w'*w)*h-w'*a
  //--------------------------------

  //help1 = w'*w
  dgemm('T', 'N', k, k, m, 1.0, *w0, m, *w0, m, 0., help1, k);
  //gradh = help1 * h
  dgemm('N', 'N', k, n, k, 1.0, help1, k, *h0, k, 0., gradh, k);
  //gradh = -w'*a + gradh
  dgemm('T', 'N', k, n, m, -1.0, *w0, m, *a, m, 1.0, gradh, k);



  double initgrad = pow(dlange('F', m, k, gradw, m, NULL), 2);
  initgrad += pow(dlange('F', k, n, gradh, k, NULL), 2);
  initgrad = sqrt(initgrad);


  double tolw = (tol > 0.001 ? tol : 0.001);
  tolw = tolw*initgrad;
  double tolh = tolw;
  


  // start of iterations

  for (iter = 1; iter <= *maxiter; ++iter) {
   

    double projnorm = 0; 
    // norm([gradW(gradW <0 | W > 0); gradH(gradH <0 | H > 0)]);
    for(i = 0; i<m*k; ++i)
      if (gradw[i] < 0 || (*w0)[i] > 0)
	projnorm += gradw[i]*gradw[i];
    for(i = 0; i<k*n; ++i)
      if (gradh[i] < 0 || (*h0)[i] > 0)
	projnorm += gradh[i]*gradh[i];
    projnorm = sqrt(projnorm);




    // stopping condition
    if (projnorm < (tol*initgrad)) {
      break;
    }


    //solve subproblems
    //-----------------



    //solve subproblem
    int iterw = pg_subprob_w(*a, *w0, *h0, gradw, hnw, hpw, help1, d, tempw, m,n,k, tolw, 1000);

    if (iterw == 1)
      tolw = 0.1 * tolw;
    
    int iterh = pg_subprob_h(*a, *w0, *h0, gradh, h, hp, help1, helph, temph, m,n,k, tolh, 1000);



    if (iterh == 1)
      tolh = 0.1 * tolh;
 



#if DEBUG_LEVEL >= 1
  dnorm = calculateNorm(*a, *w0, *h0, d, m, n, k);
  printf("iter: %.6d\t dnorm: %.16f\t proj-grad norm: %.16f\n", iter, dnorm, projnorm);
#endif   
   
  }	//end of for


    // calculating the norm of D = A-W*H
    dnorm = calculateNorm(*a, *w0, *h0, d, m, n, k);

    

    // storing the norm results of the current iteration
    dnorm0 = dnorm;



#if DEBUG_LEVEL >= 2
	printf("Exiting nmf_alspg\n");
#endif
#ifdef PROFILE_NMF_ALSPG
	gettimeofday(&end, 0);
	outputTiming("", start, end);
#endif


  // freeing memory if used
  free(gradw);
  free(gradh);
  free(help1);
  free(helph);
  free(helpw);
  free(h);
  free(hp);
  free(hn);
  free(hpw);
  free(hnw);
  free(tempw);
  free(temph);
  free(d);
  free(diffw);
  free(diffh);


  // returning calculated norm
//printf("maxiter = %d \t iter = %d \t", *maxiter, iter);
 *maxiter = iter;
  return dnorm;
}
//end of nmf_alspg
//-------------

