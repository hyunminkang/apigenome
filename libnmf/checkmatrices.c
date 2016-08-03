/**
 * checkMatrices - checks the passed matrices for non-negativity and compatibility of their dimensions
 * 
 * Purpose:
 *		This routine performs basic checks on the matrices used in nmfDriver to secure the conditions for
 *		the applied algorithms are met
 *
 * Description:
 *		This routine secures, that the matrices passed as arguments are non-negative and have compatible
 *		dimensions.
 *		matrix a - m by n
 *		matrix w - m by k
 *		matrix h - k by n
 *
 * Arguments:
 *
 * a		in, 	pointer to matrix a to check
 *
 * w		in, 	pointer to matrix w to check
 *
 * h		in, 	pointer to matrix h to check
 *
 * m		in, 	first dimension of a and w
 *
 * n		in, 	second dimension of a and h
 *
 * k		in, 	second dimension of w and first dimension of h
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <sys/time.h>
#include <time.h>

#include "common.h"
#include "outputtiming.h"


int checkMatrices(const double * a, const double * w, const double * h, const int m, \
			 const int n, const int k ) 
{

#ifdef PROFILE_CHECK_MATRICES
	struct timeval start, end;
	gettimeofday(&start, 0);
#endif
#if DEBUG_LEVEL >= 2
	printf("Entering checkMatrices\n");
#endif


  int i;
  for(i = 0; (i < m*n && a[i] >= 0.); i++) ;
  if (i < m*n) {
printf("negative element in a[%d]\n", i);
    return 1;
  }
  for(i = 0; (i < m*k && w[i] >= 0.); i++) ;
  if (i < m*k) {
printf("negative element in w[%d] = %f\n", i, w[i]);
    return 1;
  }
  for(i = 0; (i < k*n && h[i] >= 0.); i++) ;
  if (i < k*n) {
printf("negative element in h[%d]\n", i);
    return 1;
  }
#if DEBUG_LEVEL >= 2
	printf("Exiting checkMatrices\n");
#endif
#ifdef PROFILE_CHECK_MATRICES
	gettimeofday(&end, 0);
	outputTiming("Timing:", start, end);
#endif

  return 0;
}
//end of checkMatrices
//--------------------
