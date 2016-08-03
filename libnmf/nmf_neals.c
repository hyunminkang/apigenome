/**
 * NMF_NEALS -	calculates the nmf using a normal equation alternating least squares method
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
 *		The alternative least squares method as described in Berry (2006) (see references [1])
 *		is based on following algorithm (in Matlab notation):
 *		
 *		W = rand(m,k);
 *		for iter=1:maxiter
 *			H = max((W'*W)\(W'*A), 0);
 *			W = ((H*H')\(H*A'))';
 *			W = max(W, 0);
 *		end
 *
 *		The matrix right division is implemented by following Lapack routines, which solve the linear system
 *
 *		(W0'*W) * x = (W0'*A) 			(H*H' * x = A')
 *
 *		by computing a LU factorisation with partial pivoting and row interchanges
 *		
 *		A = P * L * U
 *
 *		and therefore reducing the linear system to a triangular linear system and solving
 *		that system. P is the permutation matrix of the LU factorisation.
 *		
 *		dgemm	matrix-matrix multiplications to reduce dimensionality of the problem
 *		dgesv	solving the linear system using a LU factorisation
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
#include "calculatenorm.h"
#include "calculatemaxchange.h"









double nmf_neals(double ** a, double ** w0, double ** h0, int m, int n, \
		      int k, int * maxiter, const double TolX, const double TolFun) 
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


  double * help1 = (double*) malloc(sizeof(double)*k*k);
  double * help2 = (double*) malloc(sizeof(double)*k*n);
  double * help3 = (double*) malloc(sizeof(double)*k*m);

  //-----------------------------------------



  // definition of necessary dynamic data structures
  //...for calculating matrix h
  double* h = (double*) malloc(sizeof(double)*k*n);
  int* jpvt_h = (int*) malloc(sizeof(int)*k);
  int info;
  //...for calculating matrix w
   double* w = (double*) malloc(sizeof(double)*m*k);
  //----------------


  //...for calculating the norm of A-W*H
  double* d = (double*) malloc(sizeof(double)*m*n);					//d = a - w*h
  double dnorm0 = 0;
  double dnorm = 0;
  const double eps = dlamch('E');					//machine precision epsilon
  const double sqrteps = sqrt(eps);					//squareroot of epsilon
  

  //-------------------

#ifdef ERROR_CHECKING
  if (errno) {
    perror("Error allocating memory in nmf_neals");
    free(help1);
    free(help2);
    free(help3);
    free(h);
    free(jpvt_h);
    free(w);
    free(d);
    return -1;
  }
#endif


  // declaration of data structures for switch to als algorithm
  // ----------------------------------------------------------

  int als_data_allocated = 0;					// indicates wheter data structures were already allocated
  // factor matrices for factorizing matrix w
  double * q;
  double * r;
  // factor matrices for factorizing matrix h
  double * q_h;
  double * r_h;

  double* tau_h;                   //stores elementary reflectors of factor matrix Q
  double* work_w;            //work array for factorization of matrix w
  int lwork_w;
  double* work_h;            //work array for factorization of matrix h  
  int lwork_h;
  double * work_qta;	     //work array for dorgqr
  int lwork_qta;
  double * work_qth;	     //work array for dorgqr
  int lwork_qth;

  //query for optimal workspace size for routine dgeqp3...
  double querysize;
  


  

  
  //Loop-Indices
  int iter, i;

  //variable for storing if fallback happened in current iteration
  int fallback;

  // factorisation step in a loop from 1 to maxiter
  for (iter = 1; iter <= *maxiter; ++iter) {

    //no fallback in this iteration so far
    fallback = 0;

    // calculating matrix h
    //----------------
    //help1 = w0'*w0
    dgemm('T', 'N', k, k, m, 1.0, *w0, m, *w0, m, 0., help1, k);
    //help2 = w0'*a
    dgemm('T', 'N', k, n, m, 1.0, *w0, m, *a, m, 0., help2, k);
    //LU-Factorisation of help1 to solve equation help1 * x = help2
    dgesv(k, n, help1, k, jpvt_h, help2, k, &info);
    // if factor matrix U is singular -> switch back to als algorithm to compute h
    if( info > 0) {

	//set fallback to 1 to  indicate that fallback happened
	fallback = 1;

	// do dynamic data structures need to be allocated?
	if (!als_data_allocated) {
	  als_data_allocated = 1;

  	  // factor matrices for factorizing matrix w
	  q = (double*) malloc(sizeof(double)*m*k);
	  r = (double*) malloc(sizeof(double)*m*k);
	  // factor matrices for factorizing matrix h
	  q_h = (double*) malloc(sizeof(double)*n*k);
	  r_h = (double*) malloc(sizeof(double)*n*k);

	  tau_h = (double*) malloc(sizeof(double)*k);                   //stores elementary reflectors of factor matrix Q


          //query for optimal workspace size for routine dgeqp3...
	  //for matrix w
	  dgeqp3(m, k, q, m, jpvt_h, tau_h, &querysize, -1, &info);
	  lwork_w = (int) querysize;
	  work_w = (double*) malloc(sizeof(double)*lwork_w);            //work array for factorization of matrix help1 (dgeqp3)
	  //for matrix h
	  dgeqp3(n, k, q_h, n, jpvt_h, tau_h, &querysize, -1, &info);
	  lwork_h = (int) querysize;
	  work_h = (double*) malloc(sizeof(double)*lwork_h);            //work array for factorization of matrix h


	  //query for optimal workspace size for routine dorgqr...
	  //for matrix w
	  dorgqr(m, k, k, q, m, tau_h, &querysize, -1, &info);
	  lwork_qta = (int)querysize;
	  work_qta = (double*) malloc(sizeof(double)*lwork_qta);	  //work array for dorgqr
	  //for matrix h
	  dorgqr(n, k, k, q_h, n, tau_h, &querysize, -1, &info);
	  lwork_qth = (int)querysize;
	  work_qth = (double*) malloc(sizeof(double)*lwork_qth);

	}


        // calculating matrix h
        //----------------
        //re-initialization

        //copy *w0 to q
        dlacpy('A', m, k, *w0, m, q, m);


        //initialise jpvt_h to 0 -> every column free
        for (i = 0; i<k; ++i)
          jpvt_h[i] = 0;

        // Q-R factorization with column pivoting
        dgeqp3(m, k, q, m, jpvt_h, tau_h, work_w, lwork_w, &info);


        //copying upper triangular factor-matrix r out of q into r
        dlacpy('U', m, k, q, m, r, k);


        //Begin of least-squares-solution to w0 * x = a

        //generate explicit matrix q (m times k) and calculate q' * a
        dorgqr(m, k, k, q, m, tau_h, work_qta, lwork_qta, &info);
        dgemm('T', 'N', k, n, m, 1.0, q, m, *a, m, 0.0, q_h, k);


        //solve R * x = (Q'*A)
        dtrtrs('U','N','N',k,n,r,k,q_h,k,&info);

        //copy matrix q to h, but permutated according to jpvt_h
        for (i=0; i<k; ++i) {
          dcopy(n, q_h + i, k, h + jpvt_h[i] - 1, k);
        }

        //transform negative and very small positive values to zero for performance reasons and to keep the non-negativity constraint
        for (i=0; i<k*n; ++i) {
        if (h[i] < ZERO_THRESHOLD)
          h[i] = 0.;
        }


    }
    else {
      //h = max(ZERO_THRESHOLD, help1\help2)
      for (i=0; i < k*n; ++i)
        h[i] = (help2[i] > ZERO_THRESHOLD ? help2[i] : 0.);
    }


	    // calculating matrix w = max(0, help1\help3)'
	    //----------------------------
	    //help1 = h*h'
	    dgemm('N', 'T', k, k, n, 1.0, h, k, h, k, 0., help1, k);
	    //help3 = h*a'
	    dgemm('N', 'T', k, m, n, 1.0, h, k, *a, m, 0., help3, k);
	    //LU-Factorisation of help1
	    dgesv(k, m, help1, k, jpvt_h, help3, k, &info);
	    //
	    if( info > 0) {
		// do dynamic data structures need to be allocated?
	        if (!als_data_allocated) {
	          als_data_allocated = 1;

	          // factor matrices for factorizing matrix w
	          q = (double*) malloc(sizeof(double)*m*k);
	          r = (double*) malloc(sizeof(double)*m*k);
	          // factor matrices for factorizing matrix h
	          q_h = (double*) malloc(sizeof(double)*n*k);
	          r_h = (double*) malloc(sizeof(double)*n*k);
	
	          tau_h = (double*) malloc(sizeof(double)*k);                   //stores elementary reflectors of factor matrix Q
	
	
	          //query for optimal workspace size for routine dgeqp3...
		
		  //for matrix w
        	  dgeqp3(m, k, q, m, jpvt_h, tau_h, &querysize, -1, &info);
	          lwork_w = (int) querysize;
        	  work_w = (double*) malloc(sizeof(double)*lwork_w);            //work array for factorization of matrix help1 (dgeqp3)
	
		  //..for matrix h
		  dgeqp3(n, k, q_h, n, jpvt_h, tau_h, &querysize, -1, &info);
		  lwork_h = (int) querysize;
		  work_h = (double*) malloc(sizeof(double)*lwork_h);            //work array for factorization of matrix h
	
	
        	  //query for optimal workspace size for routine dorgqr...
		  //for matrix w
	          dorgqr(m, k, k, q, m, tau_h, &querysize, -1, &info);
	          lwork_qta = (int)querysize;
	          work_qta = (double*) malloc(sizeof(double)*lwork_qta);          //work array for dorgqr
		  // ... for matrix h
		  dorgqr(n, k, k, q_h, n, tau_h, &querysize, -1, &info);
		  lwork_qth = (int)querysize;
		  work_qth = (double*) malloc(sizeof(double)*lwork_qth);
	  
	        } 

	
		//calculating matrix w
	        //copy original matrix h to q_h, but transposed
	        for (i=0; i<k; ++i) {
	          dcopy(n, h + i, k, q_h + i*n, 1);
	        }


        	//initialise jpvt_a to 0 -> every column free
	        for (i = 0; i<k; ++i)
	          jpvt_h[i] = 0;
	
        	//Q-R factorization
	        dgeqp3(n, k, q_h, n, jpvt_h, tau_h, work_h, lwork_h, &info);


	        //copying upper triangular factor-matrix r_h out of q into r_h
	        dlacpy('U', n, k, q_h, n, r_h, k);
	

	        //Begin of least-squares-solution to w0 * x = a
 
	        //generate explicit matrix q (n times k) and calculate *a = q' * a'
	        dorgqr(n, k, k, q_h, n, tau_h, work_qth, lwork_qth, &info);
	        dgemm('T', 'T', k, m, n, 1.0, q_h, n, *a, m, 0.0, q, k);




	        //solve R_h * x = (Q'*A')
	        dtrtrs('U', 'N', 'N', k, m, r_h, k, q, k, &info);


        	//jpvt_h*(R\(Q'*A')) permutation and transposed copy to w
	        for (i=0; i<k; ++i) {
	          dcopy(m, q + i, k, w + m * (jpvt_h[i] - 1), 1);
        	}

	        //transform negative and very small positive values to zero for performance reasons and to keep the non-negativity constraint
        	for (i=0; i<k*m; ++i) {
	          if (w[i] < ZERO_THRESHOLD)
	        	    w[i] = 0.;
        	}


	    }

	    else {
	        //w = max(0, help3)'
	        for (i=0; i<k; ++i) {
	          dcopy(m, help3 + i, k, w + i*m, 1);
	        }
	        for (i=0; i<m*k; ++i) {
	          if (w[i] < ZERO_THRESHOLD)
	            w[i] = 0.;
	        }
	    }




    
    // calculating the norm of D = A-W*H
    dnorm = calculateNorm(*a, w, h, d, m, n, k);

    
    // calculating change in w -> dw
    //----------------------------------
    double dw;
    dw = calculateMaxchange(w, *w0, m, k, sqrteps);

    
    // calculating change in h -> dh
    //-----------------------------------
    double dh;
    dh = calculateMaxchange(h, *h0, k, n, sqrteps);

    //Max-Change = max(dh, dw) = delta
    double delta;
    delta = (dh > dw) ? dh : dw;

   

    // storing the matrix results of the current iteration
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

    


    // storing the norm results of the current iteration
    dnorm0 = dnorm;
    


    

#if DEBUG_LEVEL >= 1
  printf("iter: %.6d\t dnorm: %.16f\t delta: %.16f\n", iter, dnorm, delta);
#endif   
   
  } //end of loop from 1 to maxiter

#if DEBUG_LEVEL >= 2
	printf("Exiting nmf_neals\n");
#endif
#ifdef PROFILE_NMF_NEALS
	gettimeofday(&end, 0);
	outputTiming("", start, end);
#endif

  // freeing memory if used
    free(help1);
    free(help2);
    free(help3);
    free(h);
    free(jpvt_h);
    free(w);
    free(d);

    if(als_data_allocated) {
      free(q);
      free(r);
      free(q_h);
      free(r_h);
      free(work_h);
      free(work_w);
      free(tau_h);
      free(work_qta);
      free(work_qth);
    }

  // returning calculated norm
  return dnorm;
}
//end of nmf_neals
//-------------
