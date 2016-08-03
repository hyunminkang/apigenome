/**
 * NMF_FAST -	calculates the nmf using a optimized normal equation least squares method
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
 *		Algorithm by Mikkel Schmidt, see Documentation for references
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





float  nmf_fast(float  ** a, float  ** w0, float  ** h0, int m, int n, \
		      int k, int * maxiter, const float  TolX, const float  TolFun) 
{
  

  
#ifdef PROFILE_NMF_FAST
	struct timeval start, end, s, e;
	gettimeofday(&start, 0);
#endif
#if DEBUG_LEVEL >= 2
	printf("Entering nmf_fast\n");
#endif


#ifdef ERROR_CHECKING
  errno = 0;
#endif

  //fastnmf data structures
  float  rowtol = 1e-8;
  int maxInnerIter = 10;
  int jc;							//current index of last unique row of aw/ah
  int row, column;
  int maxdim = (m > n) ? m : n;
  int * ic = (int*) malloc(sizeof(int) * maxdim);

  // for calculation of w
  float  * hht = (float *) malloc(sizeof(float ) * k*k);
  float  * aht = (float *) malloc(sizeof(float ) * m*k);
  float  * gw  = (float *) malloc(sizeof(float ) * m*k);   //gradient of w
  int * aw  = (int*) malloc(sizeof(int) * m*k);   	 //1s where W != 0 OR gW < 0
  int * awOrder = (int*) malloc(sizeof(int) * m); 	 // stores the order of aw rows in ascending order (interpret row as dual digit, highest bit first)
  int * awBucketSize = (int*) malloc(sizeof(int) * m);


  // for calculation of h
  float  * wtw = (float *) malloc(sizeof(float ) * k*k);
  float  * wta = (float *) malloc(sizeof(float ) * k*n);
  float  * gh  = (float *) malloc(sizeof(float ) * k*n);   //gradient of h
  int * ah  = (int*) malloc(sizeof(int) * k*n);   	 //1s where H != 0 OR gH < 0	
  int * ahOrder = (int*) malloc(sizeof(int) *n); 	 	 //stores the order of ah columns in ascending order.
  int * ahBucketSize = (int*) malloc(sizeof(int) * n);


  float  * help1 = (float *) malloc(sizeof(float )*k*k);
  float  * help2 = (float *) malloc(sizeof(float )*k*n);
  float  * help3 = (float *) malloc(sizeof(float )*k*m);
  
  int rowsHelp3, colsHelp3, rowsHelp1, colsHelp1, rowsHelp2, colsHelp2, colcount;
  //-----------------------------------------



  // definition of necessary dynamic data structures
  //...for calculating matrix h
  float * h = (float *) malloc(sizeof(float )*k*n);
  int* jpvt_h = (int*) malloc(sizeof(int)*k);
  int info;
  //...for calculating matrix w
   float * w = (float *) malloc(sizeof(float )*m*k);
  //----------------


  //...for calculating the norm of A-W*H
  float * d = (float *) malloc(sizeof(float )*m*n);					//d = a - w*h
  float  dnorm0 = 0;
  float  dnorm = 0;
  const float  eps = slamch('E');					//machine precision epsilon
  const float  sqrteps = sqrtf(eps);					//squareroot of epsilon
  

  //-------------------

#ifdef ERROR_CHECKING
  if (errno) {
    perror("Error allocating memory in nmf_fast");
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
  float  * q;
  float  * r;
  // factor matrices for factorizing matrix h
  float  * q_h;
  float  * r_h;

  float * tau_h;                   //stores elementary reflectors of factor matrix Q
  float * work_w;            //work array for factorization of matrix w
  int lwork_w;
  float * work_h;            //work array for factorization of matrix h  
  int lwork_h;
  float  * work_qta;	     //work array for sorgqr
  int lwork_qta;
  float  * work_qth;	     //work array for sorgqr
  int lwork_qth;

  //query for optimal workspace size for routine sgeqp3...
  float  querysize;
 

  
  //Loop-Indices
  int iter, i;

  //variable for storing if fallback happened in current iteration
  int fallback;
  int tind = 0;

  // factorisation step in a loop from 1 to maxiter
  for (iter = 1; iter <= *maxiter; ++iter) {

	    //no fallback in this iteration so far
	    fallback = 0;


	    // calculating matrix w = max(0, help1\help3)'
	    //----------------------------
	    //hht = h*h'
	    sgemm('N', 'T', k, k, n, 1.0, *h0, k, *h0, k, 0., hht, k);
	    //aht= a*h'
	    sgemm('N', 'T', m, k, n, 1.0, *a, m, *h0, k, 0., aht, m);
	    
	    

	    
	    //start inner loop
	    int innerIter;

	    for (innerIter = 1; innerIter <= maxInnerIter; ++innerIter) {

	      
	      //calculate gradient of W, gw = w * hht - aht;
	      slacpy('A', m, k, aht, m, gw, m);
	      sgemm('N', 'N', m, k, k, 1.0, *w0, m, hht, k, -1.0, gw, m);
	      

	      

	      float  minwgw; // to store minimum of w and gw
	      int ignorerow; // flag to check if row consists solely of entries < tol
	      int r,c, rr, cc, rowcount;		//row, column and row(uniquecheck), column(uniquecheck)
	      jc = 0;
	      //run over every row
	      for (r = 0; r < m; ++r) {
		ignorerow = 1;
		for(c = 0; c < k; ++c) {
		  minwgw = (gw[r + c * m] < (*w0)[r + c * m]) ? fabs( gw[r + c*m] ) : fabs( (*w0)[r + c*m] );

		  if (minwgw >= rowtol) {
		    ignorerow = 0;
		  }
		  if ( (*w0)[r + c*m] != 0. || gw[r + c*m] < 0.) {
		    aw[jc + c*m] = 1;				//write in next row to set
		  }
		  else {
		    aw[jc + c*m] = 0;				//write in next row to set
		  }
		}  
		//if row should be considered --> check if row jc is unique

		if( !ignorerow ) {
		  int unique = 1;
		  for(rr = 0; rr < jc; ++rr) {
		    for(cc = 0; cc < k; ++cc) {

		      if (aw[jc + cc*m] != aw[rr + cc*m])
			break;		
		    }
		    if (cc >= k) {	//non unique
		      unique = 0;
		      //store orign row number
		      ic[r] = rr;
		      break;
		    }
		  }
		  
		  //if unique --> increase jc to make current row permanent member and store origin row number
		  if (unique) {
		    //store origin row number:
		    ic[r] = jc;
		    ++jc;
		  }
		}
		else {		// if row should be ignored ->
		  ic[r] = -1;
		}
	      }
	      

	      // if no row considered --> break of innerIter loop
	      if (!jc)
		break;



	      //sort aw rows ascendingly
	      int tmpi, tmpj, tmpk, tmpl, step;
	      // start with original order awOrder[pos] = index of pos. lowest row
	      for (tmpi = 0; tmpi < jc; ++ tmpi) {
		awOrder[tmpi] = tmpi;
	      }

	      // sort aw rows in ascending order by using every column to divide the existing buckets further (starting with every row in the first bucket)
	      // final result: every row is one bucket
	      awBucketSize[0] = jc;					//index the first bucket extends to
	      for(tmpi = 1; tmpi < jc; ++tmpi)
		awBucketSize[tmpi] = 0;					//index the i-th bucket extends to from  awBucketSize[tmpi - 1]


	      // sort for every column tmpi
	      for(tmpi = 0; tmpi < k; ++tmpi) {

		//for every bucket
		int bucketstart = 0;
		for(tmpj = 0; tmpj < jc; ++tmpj) {
		  // search for first 1 in this column of this bucket and swap it with the last 0
		  int newbucketstart = awBucketSize[tmpj];
		  if  ( (awBucketSize[tmpj] - bucketstart) > 1 ) {			//if more than one row in bucket
		    for(tmpk = bucketstart; tmpk < awBucketSize[tmpj]; ++tmpk) {
		      
		      // found 1?
		      if ( aw[ awOrder[tmpk] + tmpi * m] ) {
			// start backward search from bucketend to find last 0 (to be swapped)
			for(tmpl = awBucketSize[tmpj] - 1; tmpl > tmpk; --tmpl) {
			  if ( !aw[ awOrder[tmpl] + tmpi * m]  ) {				// zero? //if ( !aw[ awOrder[tmpl] + tmpi * jc]  ) {				// zero?
			    //swap buckets and stop search for zeros.
			    int tmpPos = awOrder[tmpk];
			    awOrder[tmpk] = awOrder[tmpl];
			    awOrder[tmpl] = tmpPos;
			    
			    //stop backward search
			    break;
			  }
			} //end of backward search for last 0
			
			// no last 0 found? --> bucket is sorted regarding this column
			if (tmpl <= tmpk && bucketstart == tmpk)		// if this column cant distinguish this bucket further --> try next bucket
			  break;
			if (tmpl <= tmpk) {
			   
			  // update awBucketSize
			  // shift old bucketsizes down
			  for (tmpl = jc - 1; tmpl > tmpj; --tmpl)
			    awBucketSize[tmpl] = awBucketSize[tmpl - 1];
			  // set current Bucketsize to position of tmpk
			  awBucketSize[tmpj] = tmpk;

			  // stop search for more 1s, as bucket is now sorted in respect to current column
			  break;
			}
			
		      }
		    }
		  }
		  // update bucketstart for next iteration of loop over every bucket
		  bucketstart = newbucketstart;
		}

	      }
	    






	      // start w-update w = max(0, hht(c,c)\hat)
	      // loop over every row in aw	     
	      for(rr = 0; rr < jc; ++rr) {
		//create matrices AHt(i,c)' -> help3 and HHT(c,c) -> help1
		//use loops instead of blas/lapack routines since data is not contiguous in both dimensions
		int row, column;
		colcount = 0;
		//help3 = AHt(i,c)'
		//count rows for this submatrix = count of 1s in current aw-row
		rowsHelp3 = 0;
		for(column = 0; column < k; ++column) {
		  if ( aw[ awOrder[rr] + column * m] )
		    ++rowsHelp3;
		}
		for(row = 0; row < m; ++row) {
		  if ( ic[row] == awOrder[rr]) {
		    rowcount = 0;
		    for(column = 0; column < k; ++column) {
		      if ( aw[ awOrder[rr] + column * m ] ) {			//if ( aw[ rr + column * m] ) {
		        help3[ rowcount + colcount * rowsHelp3] = aht[ row + column * m];
			++rowcount;
		      }
		    }
		    ++colcount;

		  }
		}
  
		colsHelp3 = colcount;
		
		//help1 = HHt(c,c)
		rowcount = 0;
		for(row = 0; row < k; ++row) {
		  colcount = 0;
		  if( aw[ awOrder[rr] + row * m] ) {				//if( aw[rr + row * m] ) {
		    for(column = 0; column < k; ++column) {
		      if( aw[ awOrder[rr] + column * m] ) {				//if( aw[rr + column * m] ) {
			help1[ rowcount + colcount * rowsHelp3] = hht[row + column * k];
			++colcount;
		      }
		    }
		    ++rowcount;
		  }
		}
		
	      
		//Factorisation of hht and least square solution of system of linear equations
		//workspace query
		float  work_wupd_query = 0.;
		int lwork_wupd = 0;
		float  * work_wupd;
		ssytrf('U', rowsHelp3, help1, rowsHelp3, jpvt_h, &work_wupd_query, -1, &info);
		lwork_wupd = (int) work_wupd_query;
		work_wupd = (float *) malloc(sizeof(float ) * lwork_wupd);
		ssytrf('U', rowsHelp3, help1, rowsHelp3, jpvt_h, work_wupd, lwork_wupd, &info);
		ssytrs('U', rowsHelp3, colsHelp3, help1, rowsHelp3, jpvt_h, help3, rowsHelp3, &info);

		free(work_wupd);

		// copy result back to *w0: w0(i,c) = (help1 \ help3)' colsHelp3 x rowsHelp3
		colcount = 0;
		for(row = 0; row < m; ++row) {
		    if ( ic[row] == awOrder[rr]) {
		      rowcount = 0;
		      for(column = 0; column < k; ++column) {
			if ( aw[ awOrder[rr] + column * m] ) {
			  if ( help3[ rowcount + colcount * rowsHelp3] < ZERO_THRESHOLD)
			    (*w0)[row + column * m] = 0.;
			  else
			    (*w0)[row + column * m] = help3[ rowcount + colcount * rowsHelp3];
			  ++rowcount;
			}
		      }
		      ++colcount;
		    }
		}
	      }
	    	  
	    }


	    //
	    if( info > 0) {
		printf("!!!!!!!! fallback !!!!!!!! H-UPDATE\n");
		// do dynamic data structures need to be allocated?
	        if (!als_data_allocated) {
	          als_data_allocated = 1;

	          // factor matrices for factorizing matrix w
	          q = (float *) malloc(sizeof(float )*m*k);
	          r = (float *) malloc(sizeof(float )*m*k);
	          // factor matrices for factorizing matrix h
	          q_h = (float *) malloc(sizeof(float )*n*k);
	          r_h = (float *) malloc(sizeof(float )*n*k);
	
	          tau_h = (float *) malloc(sizeof(float )*k);                   //stores elementary reflectors of factor matrix Q
	
	
	          //query for optimal workspace size for routine sgeqp3...
		
		  //for matrix w
        	  sgeqp3(m, k, q, m, jpvt_h, tau_h, &querysize, -1, &info);
	          lwork_w = (int) querysize;
        	  work_w = (float *) malloc(sizeof(float )*lwork_w);            //work array for factorization of matrix help1 (sgeqp3)
	
		  //..for matrix h
		  sgeqp3(n, k, q_h, n, jpvt_h, tau_h, &querysize, -1, &info);
		  lwork_h = (int) querysize;
		  work_h = (float *) malloc(sizeof(float )*lwork_h);            //work array for factorization of matrix h
	
	
        	  //query for optimal workspace size for routine sorgqr...
		  //for matrix w
	          sorgqr(m, k, k, q, m, tau_h, &querysize, -1, &info);
	          lwork_qta = (int)querysize;
	          work_qta = (float *) malloc(sizeof(float )*lwork_qta);          //work array for sorgqr
		  // ... for matrix h
		  sorgqr(n, k, k, q_h, n, tau_h, &querysize, -1, &info);
		  lwork_qth = (int)querysize;
		  work_qth = (float *) malloc(sizeof(float )*lwork_qth);
	  
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


	

	    // calculating matrix h = max(0, help1\help3)'
	    //----------------------------
	    //wtw = w' * w
	    sgemm('T', 'N', k, k, m, 1.0, *w0, m, *w0, m, 0., wtw, k);
	    //wta = w' * a
	    sgemm('T', 'N', k, n, m, 1.0, *w0, m, *a, m, 0., wta, k);
	    
	    


	    //start inner loop
	    for (innerIter = 1; innerIter <= maxInnerIter; ++innerIter) {
	     
	      //calculate gradient of H, gh = wtw * h - wta;
	      slacpy('A', k, n, wta, k, gh, k);
	      sgemm('N', 'N', k, n, k, 1.0, wtw, k, *h0, k, -1.0, gh, k);

	     
      	    

	    
	     
	      

	      float  minhgh; // to store minimum of h and gh
	      int ignorecol; // flag to check if row consists solely of entries < tol
	      int r,c, rr, cc, colcount, rowcount;		//row, column and row(uniquecheck), column(uniquecheck)
	      
	      jc = 0;
	      //run over every column
	      for (c = 0; c < n; ++c) {
		ignorecol = 1;
		for(r = 0; r < k; ++r) {
		  minhgh = (gh[r + c * k] < (*h0)[r + c * k]) ? fabs( gh[r + c*k] ) : fabs( (*h0)[r + c*k] );

		  if (minhgh >= rowtol) {
		    ignorecol = 0;
		  }
		  if ( (*h0)[r + c*k] != 0. || gh[r + c*k] < 0.)
		    ah[r + jc * k] = 1;				//write in row next to be set
		  else
		    ah[r + jc * k] = 0;				//write in row next to be set
		}  
		//if column should be considered --> check if column jc is unique
		
		if( !ignorecol ) {
		  int unique = 1;
		  for(cc = 0; cc < jc; ++cc) {
		    for(rr = 0; rr < k; ++rr) {
		      if (ah[rr + jc * k] != ah[rr + cc*k])
			break;		
		    }
		    if (rr >= k) {	//non unique 
		      unique = 0;
		      //store orign col number
		      ic[c] = cc;
		      break;
		    }
		  }
		  
		  //if unique --> increase jc to make current row permanent member and store origin row number
		  if (unique) {
		    //store origin row number:
		    ic[c] = jc;
		    ++jc;
		  }
		}
		else { 		// if col should be ignored
		  ic[c] = -1;
		}  
	      }
	      
	      

	      // if no row considered --> break of innerIter loop
	      if (!jc)
		break;

	      
	      
	      //sort ah columns ascendingly
	      int tmpi, tmpj, tmpk, tmpl, step;
	      // start with original order ahOrder[pos] = index of pos. lowest column
	      for (tmpi = 0; tmpi < jc; ++ tmpi) {
		ahOrder[tmpi] = tmpi;
	      }

	      // sort ah cols in ascending order by using every row to divide the existing buckets further (starting with every col in the first bucket)
	      // final result: every col is one bucket

	      ahBucketSize[0] = jc;					//index the first bucket extends to
	      for(tmpi = 1; tmpi < jc; ++tmpi)
		ahBucketSize[tmpi] = 0;					//index the i-th bucket extends to from  awBucketSize[tmpi - 1]

	      // sort for every row tmpi
	      for(tmpi = 0; tmpi < k; ++tmpi) {
		
		//for every bucket
		int bucketstart = 0;
		for(tmpj = 0; tmpj < jc; ++tmpj) {
		  // search for first 1 in this column of this bucket and swap it with the last 0
		  int newbucketstart = awBucketSize[tmpj];
		  if ( (awBucketSize[tmpj] - bucketstart) > 1) {
		    for(tmpk = bucketstart; tmpk < ahBucketSize[tmpj]; ++tmpk) {
		      
		      // found 1?
		      if ( ah[ tmpi + ahOrder[tmpk] * k] ) {
			// start backward search from bucketend to find last 0 (to be swapped)
			for(tmpl = ahBucketSize[tmpj] - 1; tmpl > tmpk; --tmpl) {
			  if ( !ah[ tmpi + ahOrder[tmpl] * k]  ) {				// zero?
			    //swap buckets and stop search for zeros.
			    int tmpPos = ahOrder[tmpk];
			    ahOrder[tmpk] = ahOrder[tmpl];
			    ahOrder[tmpl] = tmpPos;
			    
			    //stop backward search
			    break;
			  }
			} //end of backward search for last 0
			
			// no last 0 found? --> bucket is sorted regarding this column 
			if (tmpl <= tmpk && bucketstart == tmpk)		// if this column cant distinguish this bucket further --> try next bucket
			  break;
			if (tmpl <= tmpk) {
			    
			  // update awBucketSize
			  // shift old bucketsizes down
			  for (tmpl = jc - 1; tmpl > tmpj; --tmpl)
			    ahBucketSize[tmpl] = ahBucketSize[tmpl - 1];
			  // set current Bucketsize to position of tmpk
			  ahBucketSize[tmpj] = tmpk;
			  // stop search for more 1s, as bucket is now sorted in respect to current column
			  break;
			}
			
		      }
		    }// end of for tmpk....
		  }//enf of if 
		  // update bucketstart for next iteration of loop over every bucket
		  bucketstart = newbucketstart;
		}
	      }	      
	      


	      
	      
	      

	      // start h-update h = max(0, wtw \ wta)
	      // loop over every col in ah
	      for(cc = 0; cc < jc; ++cc) {
		//create matrices wtw(c,c) -> help1 and wta(c,j) -> help2
		//use loops instead of blas/lapack routines since data is not contiguous in both dimensions
		int row, column;
		rowcount = 0;
		//help2 = wta(c,j)
		//count rows for this submatrix = count of 1s in current ah-column
		rowsHelp2 = 0;
		for(row = 0; row < k; ++row) {
		  if ( ah[row + ahOrder[cc] * k] )
		    ++rowsHelp2;
		}
		for(row = 0; row < k; ++row) {
		  if ( ah[ row + ahOrder[cc] * k]) { 
		    colcount = 0;
		    for(column = 0; column < n; ++column) {
		      if ( ic[column] == ahOrder[cc]) {
		        help2[ rowcount + colcount * rowsHelp2] = wta[row + column * k];
			++colcount;
		      }
		    }
		    ++rowcount;
		  }
		}
		colsHelp2 = colcount;
		
	
		//help1 = WtW(c,c)
		rowcount = 0;
		for(row = 0; row < k; ++row) {
		  colcount = 0;
		  if( ah[row + ahOrder[cc] * k] ) {
		    for(column = 0; column < k; ++column) {
		      if( ah[column + ahOrder[cc] * k] ) {
			help1[ rowcount + colcount * rowsHelp2] = wtw[row + column * k];
			++colcount;
		      }
		    }
		    ++rowcount;
		  }
		}
		
		//Factorisation of wtw
		float  work_hupd_query = 0.;
		int lwork_hupd = 0;
		float  *work_hupd;
		ssytrf('U', rowsHelp2, help1, rowsHelp2, jpvt_h, &work_hupd_query, -1, &info);
		lwork_hupd = (int) work_hupd_query;
		work_hupd = (float *) malloc(sizeof(float ) * lwork_hupd);
		ssytrf('U', rowsHelp2, help1, rowsHelp2, jpvt_h, work_hupd, lwork_hupd, &info);
		ssytrs('U', rowsHelp2, colsHelp2, help1, rowsHelp2, jpvt_h, help2, rowsHelp2, &info);

		free(work_hupd);
		

		// copy result back to *h0: h0(i,c) = (help1 \ help2) rowsHelp2 x colsHelp2
		rowcount = 0;
		for(row = 0; row < k; ++row) {
		    if ( ah[row + ahOrder[cc] * k] ) {
		      colcount = 0;
		      for(column = 0; column < n; ++column) {
			if ( ic[column] == ahOrder[cc] ) {
			  if ( help2[ rowcount + colcount * rowsHelp2] < ZERO_THRESHOLD)
			    (*h0)[row + column * k] = 0.;
			  else
			    (*h0)[row + column * k] = help2[ rowcount + colcount * rowsHelp2];
			  ++colcount;
			}
		      }
		      ++rowcount;
		    }
		}
	      }

	    }
	    


    if( info > 0) {
	printf("!!!!!!!!! fallback !!!!!!!!!!!! W-UPDATE");
	//set fallback to 1 to  indicate that fallback happened
	fallback = 1;

	// do dynamic data structures need to be allocated?
	if (!als_data_allocated) {
	  als_data_allocated = 1;

  	  // factor matrices for factorizing matrix w
	  q = (float *) malloc(sizeof(float )*m*k);
	  r = (float *) malloc(sizeof(float )*m*k);
	  // factor matrices for factorizing matrix h
	  q_h = (float *) malloc(sizeof(float )*n*k);
	  r_h = (float *) malloc(sizeof(float )*n*k);

	  tau_h = (float *) malloc(sizeof(float )*k);                   //stores elementary reflectors of factor matrix Q


          //query for optimal workspace size for routine sgeqp3...
	  //for matrix w
	  sgeqp3(m, k, q, m, jpvt_h, tau_h, &querysize, -1, &info);
	  lwork_w = (int) querysize;
	  work_w = (float *) malloc(sizeof(float )*lwork_w);            //work array for factorization of matrix help1 (sgeqp3)
	  //for matrix h
	  sgeqp3(n, k, q_h, n, jpvt_h, tau_h, &querysize, -1, &info);
	  lwork_h = (int) querysize;
	  work_h = (float *) malloc(sizeof(float )*lwork_h);            //work array for factorization of matrix h


	  //query for optimal workspace size for routine sorgqr...
	  //for matrix w
	  sorgqr(m, k, k, q, m, tau_h, &querysize, -1, &info);
	  lwork_qta = (int)querysize;
	  work_qta = (float *) malloc(sizeof(float )*lwork_qta);	  //work array for sorgqr
	  //for matrix h
	  sorgqr(n, k, k, q_h, n, tau_h, &querysize, -1, &info);
	  lwork_qth = (int)querysize;
	  work_qth = (float *) malloc(sizeof(float )*lwork_qth);

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
    
    


  

#if DEBUG_LEVEL >= 1
  printf("iter: %.6d\t dnorm: %.16f\t delta: %.16f\n", iter, dnorm, delta);
#endif   
   
  } //end of loop from 1 to maxiter


    // calculating the norm of D = A-W*H
    dnorm = calculateNorm_singleprec(*a, *w0, *h0, d, m, n, k);
    
#if DEBUG_LEVEL >= 2
	printf("Exiting nmf_fast\n");
#endif
#ifdef PROFILE_NMF_FAST
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
    free(awBucketSize);
    free(aw);
    free(awOrder);
    free(ah);
    free(ahOrder);
    free(ahBucketSize);
    free(wta);
    free(wtw);
    free(hht);
    free(aht);
    free(gh);
    free(gw);
    free(ic);

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
