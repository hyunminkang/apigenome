/**
 * calculateMaxchange - calculates MaxChange of matrices w or h in the current iteration
 *
 * Purpose:
 *		Calculates the maximum change matrix w or h show from last iteration to the current iteration, which
 *		is used as one of the convergence checks in most implemented nmf algorithms
 *
 * Description:	
 *		Maximum change is calculated as max(max(abs(mat-mat0) / (sqrteps+max(max(abs(mat0))))))
 *		Uses Lapack routine "dlange" to compute the absolute maximum element of a matrix
 *
 * Arguments:
 *
 * mat			in, 	pointer to matrix of the current iteration
 *
 * mat0			in, 	pointer to matrix of the last iteration
 *				on exit mat0 = mat0 - mat
 *
 * m			in, 	first dimension of matrices mat, mat0
 *
 * n			in, 	second dimension of matrices mat, mat0
 *
 * sqrteps		in, 	sqareroot of the machine precision epsilon
 *
 * function value	out, 	calculated maxChange
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



double calculateMaxchange(double * mat, double * mat0, int m, int n, const double sqrteps)
{

#ifdef PROFILE_CALCULATE_MAXCHANGE
	struct timeval start, end;
	gettimeofday(&start, 0);
#endif
#if DEBUG_LEVEL >= 2
	printf("Entering calculateMaxchange\n");
#endif


    //absolute maximum element in m0
    double temp = dlange('M', m, n, mat0, m, NULL);

    //calculating mat0 = mat0-mat
    daxpy(m*n, - 1.0, mat, 1, mat0, 1);
    //dmat = max(max(abs(mat-mat0) / (sqrteps+max(max(abs(mat0))))));
    double dmat = dlange('M', m, n, mat0, m, NULL) / (sqrteps + temp);		

#if DEBUG_LEVEL >= 2
	printf("Exiting calculateMaxchange\n");
#endif
#ifdef PROFILE_CALCULATE_MAXCHANGE
	gettimeofday(&end, 0);
	outputTiming("Timing:", start, end);
#endif

    return dmat;
}
//end of calculateMaxchange
//-------------------------
