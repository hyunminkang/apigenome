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
#include "calculatenorm_singleprec.h"
#include "calculatemaxchange_singleprec.h"









float nmf_neals_singleprec(float ** a, float ** w0, float ** h0, int m, int n, \
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

  //-----------------------------------------



  // definition of necessary dynamic data structures
  //...for calculating matrix h
  float* h = (float*) malloc(sizeof(float)*k*n);
  int* jpvt_h = (int*) malloc(sizeof(int)*k);
  int info;
  //...for calculating matrix w
   float* w = (float*) malloc(sizeof(float)*m*k);
  //----------------


  //...for calculating the norm of A-W*H
  float* d = (float*) malloc(sizeof(float)*m*n);					//d = a - w*h
  float dnorm0 = 0;
  float dnorm = 0;
  const float eps = slamch('E');					//machine precision epsilon
  const float sqrteps = sqrt(eps);					//squareroot of epsilon
  

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
  float * q;
  float * r;
  // factor matrices for factorizing matrix h
  float * q_h;
  float * r_h;

  float* tau_h;                   //stores elementary reflectors of factor matrix Q
  float* work_w;            //work array for factorization of matrix w
  int lwork_w;
  float* work_h;            //work array for factorization of matrix h  
  int lwork_h;
  float * work_qta;	     //work array for dorgqr
  int lwork_qta;
  float * work_qth;	     //work array for dorgqr
  int lwork_qth;

  //query for optimal workspace size for routine dgeqp3...
  float querysize;
  


  

  
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
    sgemm('T', 'N', k, k, m, 1.0, *w0, m, *w0, m, 0., help1, k);
    //help2 = w0'*a
    sgemm('T', 'N', k, n, m, 1.0, *w0, m, *a, m, 0., help2, k);
    //LU-Factorisation of help1 to solve equation help1 * x = help2
    sgesv(k, n, help1, k, jpvt_h, help2, k, &info);
    // if factor matrix U is singular -> switch back to als algorithm to compute h
    if( info > 0) {

	//set fallback to 1 to  indicate that fallback happened
	fallback = 1;

	// do dynamic data structures need to be allocated?
	if (!als_data_allocated) {
	  als_data_allocated = 1;

  	  // factor matrices for factorizing matrix w
	  q = (float*) malloc(sizeof(float)*m*k);
	  r = (float*) malloc(sizeof(float)*m*k);
	  // factor matrices for factorizing matrix h
	  q_h = (float*) malloc(sizeof(float)*n*k);
	  r_h = (float*) malloc(sizeof(float)*n*k);

	  tau_h = (float*) malloc(sizeof(float)*k);                   //stores elementary reflectors of factor matrix Q


          //query for optimal workspace size for routine dgeqp3...
	  //for matrix w
	  sgeqp3(m, k, q, m, jpvt_h, tau_h, &querysize, -1, &info);
	  lwork_w = (int) querysize;
	  work_w = (float*) malloc(sizeof(float)*lwork_w);            //work array for factorization of matrix help1 (dgeqp3)
	  //for matrix h
	  sgeqp3(n, k, q_h, n, jpvt_h, tau_h, &querysize, -1, &info);
	  lwork_h = (int) querysize;
	  work_h = (float*) malloc(sizeof(float)*lwork_h);            //work array for factorization of matrix h


	  //query for optimal workspace size for routine dorgqr...
	  //for matrix w
	  sorgqr(m, k, k, q, m, tau_h, &querysize, -1, &info);
	  lwork_qta = (int)querysize;
	  work_qta = (float*) malloc(sizeof(float)*lwork_qta);	  //work array for dorgqr
	  //for matrix h
	  sorgqr(n, k, k, q_h, n, tau_h, &querysize, -1, &info);
	  lwork_qth = (int)querysize;
	  work_qth = (float*) malloc(sizeof(float)*lwork_qth);

	}


        // calculating matrix h
        //----------------
        //re-initialization

        //copy *w0 to q
        slacpy('A', m, k, *w0, m, q, m);


        //initialise jpvt_h to 0 -> every column free
        for (i = 0; i<k; ++i)
          jpvt_h[i] = 0;

        // Q-R factorization with column pivoting
        sgeqp3(m, k, q, m, jpvt_h, tau_h, work_w, lwork_w, &info);


        //copying upper triangular factor-matrix r out of q into r
        slacpy('U', m, k, q, m, r, k);


        //Begin of least-squares-solution to w0 * x = a

        //generate explicit matrix q (m times k) and calculate q' * a
        sorgqr(m, k, k, q, m, tau_h, work_qta, lwork_qta, &info);
        sgemm('T', 'N', k, n, m, 1.0, q, m, *a, m, 0.0, q_h, k);


        //solve R * x = (Q'*A)
        strtrs('U','N','N',k,n,r,k,q_h,k,&info);

        //copy matrix q to h, but permutated according to jpvt_h
        for (i=0; i<k; ++i) {
          scopy(n, q_h + i, k, h + jpvt_h[i] - 1, k);
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
	    sgemm('N', 'T', k, k, n, 1.0, h, k, h, k, 0., help1, k);
	    //help3 = h*a'
	    sgemm('N', 'T', k, m, n, 1.0, h, k, *a, m, 0., help3, k);
	    //LU-Factorisation of help1
	    sgesv(k, m, help1, k, jpvt_h, help3, k, &info);
	    //
	    if( info > 0) {
		// do dynamic data structures need to be allocated?
	        if (!als_data_allocated) {
	          als_data_allocated = 1;

	          // factor matrices for factorizing matrix w
	          q = (float*) malloc(sizeof(float)*m*k);
	          r = (float*) malloc(sizeof(float)*m*k);
	          // factor matrices for factorizing matrix h
	          q_h = (float*) malloc(sizeof(float)*n*k);
	          r_h = (float*) malloc(sizeof(float)*n*k);
	
	          tau_h = (float*) malloc(sizeof(float)*k);                   //stores elementary reflectors of factor matrix Q
	
	
	          //query for optimal workspace size for routine dgeqp3...
		
		  //for matrix w
        	  sgeqp3(m, k, q, m, jpvt_h, tau_h, &querysize, -1, &info);
	          lwork_w = (int) querysize;
        	  work_w = (float*) malloc(sizeof(float)*lwork_w);            //work array for factorization of matrix help1 (dgeqp3)
	
		  //..for matrix h
		  sgeqp3(n, k, q_h, n, jpvt_h, tau_h, &querysize, -1, &info);
		  lwork_h = (int) querysize;
		  work_h = (float*) malloc(sizeof(float)*lwork_h);            //work array for factorization of matrix h
	
	
        	  //query for optimal workspace size for routine dorgqr...
		  //for matrix w
	          sorgqr(m, k, k, q, m, tau_h, &querysize, -1, &info);
	          lwork_qta = (int)querysize;
	          work_qta = (float*) malloc(sizeof(float)*lwork_qta);          //work array for dorgqr
		  // ... for matrix h
		  sorgqr(n, k, k, q_h, n, tau_h, &querysize, -1, &info);
		  lwork_qth = (int)querysize;
		  work_qth = (float*) malloc(sizeof(float)*lwork_qth);
	  
	        } 

	
		//calculating matrix w
	        //copy original matrix h to q_h, but transposed
	        for (i=0; i<k; ++i) {
	          scopy(n, h + i, k, q_h + i*n, 1);
	        }


        	//initialise jpvt_a to 0 -> every column free
	        for (i = 0; i<k; ++i)
	          jpvt_h[i] = 0;
	
        	//Q-R factorization
	        sgeqp3(n, k, q_h, n, jpvt_h, tau_h, work_h, lwork_h, &info);


	        //copying upper triangular factor-matrix r_h out of q into r_h
	        slacpy('U', n, k, q_h, n, r_h, k);
	

	        //Begin of least-squares-solution to w0 * x = a
 
	        //generate explicit matrix q (n times k) and calculate *a = q' * a'
	        sorgqr(n, k, k, q_h, n, tau_h, work_qth, lwork_qth, &info);
	        sgemm('T', 'T', k, m, n, 1.0, q_h, n, *a, m, 0.0, q, k);




	        //solve R_h * x = (Q'*A')
	        strtrs('U', 'N', 'N', k, m, r_h, k, q, k, &info);


        	//jpvt_h*(R\(Q'*A')) permutation and transposed copy to w
	        for (i=0; i<k; ++i) {
	          scopy(m, q + i, k, w + m * (jpvt_h[i] - 1), 1);
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
	          scopy(m, help3 + i, k, w + i*m, 1);
	        }
	        for (i=0; i<m*k; ++i) {
	          if (w[i] < ZERO_THRESHOLD)
	            w[i] = 0.;
	        }
	    }




    
    // calculating the norm of D = A-W*H
    dnorm = calculateNorm_singleprec(*a, w, h, d, m, n, k);

    
    // calculating change in w -> dw
    //----------------------------------
    float dw;
    dw = calculateMaxchange_singleprec(w, *w0, m, k, sqrteps);

    
    // calculating change in h -> dh
    //-----------------------------------
    float dh;
    dh = calculateMaxchange_singleprec(h, *h0, k, n, sqrteps);

    //Max-Change = max(dh, dw) = delta
    float delta;
    delta = (dh > dw) ? dh : dw;

   

    // storing the matrix results of the current iteration
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
