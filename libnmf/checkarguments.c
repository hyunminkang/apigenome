/**
 * checkArguments - performs basic checks on the arguments passed to nmfDriver
 *
 * Purpose:
 *		Performs sanity checks on the arguments passed to nmfDriver to check for inconsistencies and errors
 *		of the passed arguments
 *
 * Arguments:
 *
 * a		in, 	filename to load matrix a from
 *
 * k		in, 	approximation factor
 *
 * iter		in, 	maximal number of iterations to perform
 *
 * rep		in, 	number of repeated factorisations with differently initialised matrices w0 and h0
 *
 * alg		in, 	algorithm to use for the factorisation
 *
 * init		in, 	defines how to initialise the matrices w0 and h0
 *
 * min_init	in, 	defines the minimum value for initialisation
 *
 * max_init	in, 	defines the maximum value for initialisation
 *
 * w0		in, 	filename to load matrix w0 from
 *
 * h0		in, 	filename to load matrix h0 from
 *
 * w_out	in, 	filename to store final matrix w to
 *
 * h_out	in, 	filename to store final matrix h to
 *
 * TolX		in, 	tolerance value for convergence check of maxChange
 *
 * TolFun	in, 	tolerance value for convergence check of dnorm
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>

#include "common.h"
#include "outputtiming.h"


int checkArguments(const char* a, const int k, int iter, const char* w0, const char* h0, options_t * opts)
{

#ifdef PROFILE_CHECK_ARGUMENTS
	struct timeval start, end;
	gettimeofday(&start, 0);
#endif
#if DEBUG_LEVEL >= 2
	printf("Entering checkArguments\n");
#endif

  if (!a || k < 0 || iter < 0 || (w0 && w0[0] == '\0' ) || (h0 && h0[0] == '\0') || opts->rep <0 || opts->min_init < 0 || !opts->w_out || !opts->h_out || opts->TolX < 0 || opts->TolFun < 0) {
    errno = EDOM;
    perror("Error in arguments passed to nmfDriver.");
    return 1;
  }



#if DEBUG_LEVEL >= 2
	printf("Exiting checkArguments\n");
#endif
#ifdef PROFILE_CHECK_ARGUMENTS
	gettimeofday(&end, 0);
	outputTiming("Timing:", start, end);
#endif

  return 0;
}
//end of checkArguments
//---------------------
