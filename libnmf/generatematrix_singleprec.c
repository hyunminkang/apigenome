/**
 * generateMatrix - generates initialization matrices W0 and H0
 *
 * Purpose
 *		Allocating memory for both matrices and initialising them using the method specified by "init"
 *
 * Description
 *		Supports following initialization strategies:
 *			Random		random numbers, min and max can be set in options_t structure "opts" passed to the routine
 *					default-option, if "null" is passed as parameter "opts" (range 0 - 1)
 *					@see options_t in common.h
 * 					enum value "init_rand"
 *
 *			NNDSVD		non negative float  singular value decomposition
 *					option can be set using options_t structure "opts" passed to the routine
 *					@see options_t in common.h
 * 					enum value "init_nndsvd"
 * 
 *			InfoGain	Feature selection based on information gain of all attributes
 *					vector of class labels and count of distinct class values needs to be passed in options_t structure
 *					@see options_t in common.h 
 * 					enum value "init_infogain"
 *
 *			GainRatio	Feature selection based on gain ratio of all attributes
 *					vector of class labels and count of distinct class values needs to be passed in options_t structure
 *					@see options_t in common.h 
 * 					enum value "init_gainratio"
 * 
 * Arguments:
 *
 * m		in, 	first dimension of matrix
 *
 * n		in, 	second dimension of matrix 
 *
 * k		in,	approximation factor
 *
 * init		in, 	type of initialization
 *
 * min		in, 	lower bound of random numbers
 *
 * max		in, 	upper bound of random numbers
 *
 * matrixW	in/out, pointer to memory storing the matrix (m x k)
 *
 * matrixH	in/out, pointer to memory storing the matrix (k x n)
 * 
 * matrixA	in,	pointer to matrix which should be factorised
 *
 * opts		in,	pointer to options_t structure to set initialization options, @see "options_t" in common.h
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
#include "randnumber_singleprec.h"
#include "calculatesvd_singleprec.h"
#include "infogain_singleprec.h"
#include "infogain_openmp_singleprec.h"





void generateMatrix_singleprec(const int m, const int n, const int k, init_t init, const int min, const int max, float  **matrixW, float  **matrixH, float  * matrixA, options_t_singleprec * opts)
{

#ifdef PROFILE_GENERATE_MATRIX
	struct timeval start, end;
	gettimeofday(&start, 0);
#endif
#if DEBUG_LEVEL >= 2
	printf("Entering generateMatrix\n");
#endif


  //if *matrix is NULL, memory has to be allocated
  if ( !(*matrixW) )
    *matrixW = (float *) malloc(sizeof(float )*m*k);
  if ( !(*matrixH) )      
    *matrixH = (float *) malloc(sizeof(float )*k*n);

#ifdef ERROR_CHECKING
  if (errno) {
    perror("Failed to allocate memory in generateMatrix");
    return;
  }
#endif

  if (init == init_infogain || init ==  init_gainratio) {
    //check for requirements
    if (opts->num_dist_classes < 1 || opts->classColumn == NULL) {
	printf("Warning: illegal classColumn and count of distinct classes, falling back to random initialization\n");
	init = init_rand;
    }
    else {
      int i;					//loop variable
      //calculate infogain and gainratio
      idx_float  * indexedInfoGainSorted = (idx_float *) malloc(sizeof(idx_float ) * n);
      idx_float  * indexedGainRatioSorted = (idx_float *) malloc(sizeof(idx_float ) * n);
    
      int error;
      if (opts->parallelization == par_openmp) {
	error = infoGain_openmp_singleprec(matrixA, m, n, opts->classColumn, opts->num_dist_classes, indexedInfoGainSorted, indexedGainRatioSorted);
      }
      else {
	error = infoGain_singleprec(matrixA, m, n, opts->classColumn, opts->num_dist_classes, indexedInfoGainSorted, indexedGainRatioSorted);
      }
     
      //copy k first (wrt information gain or gain ratio) columns of A to matrixW
      if (init == init_infogain) {
	for (i = 0; i < k; ++i) {
	  int srcRow = indexedInfoGainSorted[i].idx;
	  scopy(m, matrixA + srcRow * m, 1, (*matrixW) + i * m, 1);
	}
      }
      else {
	for (i = 0; i < k; ++i) {
	  int srcRow = indexedGainRatioSorted[i].idx;
	  scopy(m, matrixA + srcRow * m, 1, (*matrixW) + i * m, 1);
	}	
      }
      
      
      //free memory
      free(indexedInfoGainSorted);
      free(indexedGainRatioSorted);
      
      //randomly initialize factor matrix H
      for(i = 0; i < k*n; ++i)
	(*matrixH)[i] = randnumber_singleprec(min, max);
    }
  }  
  //random initialization
  if (init == init_rand) {
    int i;
    for(i = 0; i < m*k; ++i)
      (*matrixW)[i] = randnumber_singleprec(min, max);
    for(i = 0; i < k*n; ++i)
      (*matrixH)[i] = randnumber_singleprec(min, max);
  }
  if (init == init_nndsvd) {
    float  * U, *S, *V;
    
    int mindim = (m > n) ? n : m;
    int vcols, ucols;
    
    if (opts->nndsvd_maxiter < 0)
      opts->nndsvd_maxiter = 300;
    if (opts->nndsvd_ncv < 0) {
      if ( k > 10) {
        opts->nndsvd_ncv = 2 * k;
      }
      else {
	opts->nndsvd_ncv = 20;
      }
    }
    if (opts->nndsvd_ncv > mindim)
      opts->nndsvd_ncv = mindim;
    if (opts->nndsvd_ncv <= k)
      opts->nndsvd_ncv = k + 1;
   

    if (m >= n) {
      vcols = opts->nndsvd_ncv;
      ucols = k;
    }
    else {
      vcols = k;
      ucols = opts->nndsvd_ncv;
    }
    float  *uup, *uun, *vvp, *vvn;
    U = (float *) malloc(sizeof(float ) * m * ucols);
    S = (float *) malloc(sizeof(float ) * k);
    V = (float *) malloc(sizeof(float ) * n * vcols);
    uup = (float *) malloc(sizeof(float ) * m);
    uun = (float *) malloc(sizeof(float ) * m);
    vvp = (float *) malloc(sizeof(float ) * n);
    vvn = (float *) malloc(sizeof(float ) * n);
    float  n_uup, n_uun, n_vvp, n_vvn, termp, termn;
    opts->num_dist_classes = -1;
    opts->classColumn = NULL;
    

    //calculate svd
    int errcode = 0;
    errcode = calculateSVD_singleprec(matrixA, U, S, V, m, n, k, opts->nndsvd_maxiter, opts->nndsvd_tol, opts->nndsvd_ncv, 1);

#ifdef ERROR_CHECKING

    if (errcode < 0) {
	    errno = errcode;
	    free(U);
	    free(S);
	    free(V);
	    free(uup);
	    free(uun);
	    free(vvp);
	    free(vvn);
	    perror("Error in calculateSVD. Aborting.");
	    return;
    }
#endif

    if (errcode < k)
	printf("Warning: Only %d singular values of %d calculated converged\n", errcode, k);
    

    
    //fill outputmatrices with zeros
    slaset('A', m, k, 0., 0., *matrixW, m);
    slaset('A', k, n, 0., 0., *matrixH, k);
    
    //setting first column of W and first row of H
    float  sqrt_sv = sqrt(S[0]);
    saxpy(m, sqrt_sv, U, 1, *matrixW, 1);
    saxpy(n, sqrt_sv, V, 1, *matrixH, k);

    
    
    //setting remaining columns of W and rows of H
    int i, j;
    for(i = 1; i < k; ++i) {

	for (j = 0; j < m; ++j) {
		if (U[j + i*m] < 0.) {
			uup[j] = 0.0;
			uun[j] = fabsf(U[j + i*m]);
		}
		else {
			uun[j] = 0.0;
			uup[j] = U[j + i*m];
		}
	}

	for (j = 0; j < n; ++j) {
		if (V[j + i*n] < 0.) {
			vvp[j] = 0.0;
			vvn[j] = fabsf(V[j + i*n]);
		}
		else {
			vvn[j] = 0.0;
			vvp[j] = V[j + i*n];
		}
	}

	n_uup = snrm2(m, uup, 1);
	n_uun = snrm2(m, uun, 1);
	n_vvp = snrm2(n, vvp, 1);
	n_vvn = snrm2(n, vvn, 1);
	termp = n_uup * n_vvp;
	termn = n_uun * n_vvn;

	if (termp >= termn) {
		float  factorw = sqrt(S[i] * termp) / n_uup;
		float  factorh = sqrt(S[i] * termp) / n_vvp;
		saxpy(m, factorw, uup, 1, *matrixW + i*m, 1);
		saxpy(n, factorh, vvp, 1, *matrixH + i, k); 
	}
	else {
		float  factorw = sqrt(S[i] * termn) / n_uun;
		float  factorh = sqrt(S[i] * termn) / n_vvn;
		saxpy(m, factorw, uun, 1, *matrixW + i*m, 1);
		saxpy(n, factorh, vvn, 1, *matrixH + i, k);

	}
	
      
    }
    
    //set negative elements to absolute value or zero if absolute value is smaller than ZERO_THRESHOLD from common.h
    for(i = 0; i < m*k; ++i) {
      float  tmp_abs = 0.0;
      tmp_abs = fabsf((*matrixW)[i]);
      if (tmp_abs < ZERO_THRESHOLD) {
 	tmp_abs = 0.0;
      }
      (*matrixW)[i] = tmp_abs;
    }
    
    //set negative elements to absolute value or zero if absolute value is smaller than ZERO_THRESHOLD from common.h
    for(i = 0; i < k*n ; ++i) {
      float  tmp_abs = 0.0;
      tmp_abs = fabsf((*matrixH)[i]);
      if (tmp_abs < ZERO_THRESHOLD) {
 	tmp_abs = 0.0;
      }
      (*matrixH)[i] = tmp_abs;    
    }

    
    free(U);
    free(S);
    free(V);
    free(uup);
    free(uun);
    free(vvp);
    free(vvn);
  }


#if DEBUG_LEVEL >= 2
	printf("Exiting generateMatrix\n");
#endif
#ifdef PROFILE_GENERATE_MATRIX
	gettimeofday(&end, 0);
	outputTiming("", start, end);
#endif
  return;
}
//end of generateMatrix
//---------------------
