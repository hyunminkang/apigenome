/**
 * calculateNorm_singleprec - calculates the norm of A-W*H in single precision
 *
 * Purpose:
 *		Calculates the frobenius norm of matrix d = A-W*H
 *		The frobenius norm is then divided by sqrt(m*n) for normalization
 *
 * Description:
 *		This function uses the Lapack routine "slange"
 *
 * Arguments:
 *
 * a			in, 	pointer to Matrix a
 *
 * w			in, 	pointer to Matrix w
 *
 * h			in, 	pointer to Matrix h
 *
 * d			in/out, pointer to allocated storage space for matrix d=a-w*h
 *				on exit d = a - w*h
 * m			in, 	first dimension of matrix a and w
 *
 * n			in, 	second dimension of matrix a and h
 *
 * k			in, 	approximation factor, second dimension of matrix w, first dimension of matrix h
 *
 * function value	out, 	frobenius norm of A-W*H / sqrt(m*n)
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



float calculateNorm_singleprec(float *a, float *w, float *h, float * d, int m, int n, int k)
{

#ifdef PROFILE_CALCULATE_NORM
	struct timeval start, end;
	gettimeofday(&start, 0);
#endif
#if DEBUG_LEVEL >= 2
	printf("Entering calculateNorm\n");
#endif

  
  float dnorm = 0;

  slacpy('A', m, n, a, m, d, m);

  //d = - w*h
  sgemm('N', 'N', m, n, k, - 1.0, w, m, h, k, 1., d, m);


  //dnorm = frobnorm(d) / sqrt(m*n)
  dnorm = slange('F', m, n, d, m, NULL);
  dnorm = dnorm / sqrt(m*n);
  

#if DEBUG_LEVEL >= 2
	printf("Exiting calculateNorm\n");
#endif
#ifdef PROFILE_CALCULATE_NORM
	gettimeofday(&end, 0);
	outputTiming("Timing:", start, end);
#endif

  return dnorm;
}
//end of calculateNorm
//--------------------
