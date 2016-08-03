/** storeMatrix
* Purpose:
 * 		Stores a matrix to the specified file
 *
 * Description:
 * 		The file has to be a simple ascii file
 *		The first row has exactly one entry which is the number of rows m
 *		The second row has exactly one entry which is the number of columns n
 *		The next m rows contain the elements of the matrix in row-major-order
 *		Every element is separated by a single space
 *
 * Arguments:
 *
 * fileName	in, 	The path to the file to write
 *
 * m		in, 	first dimension of the stored matrix
 *
 * n		in, 	second dimension of the stored matrix
 *
 * matrix	in, 	m*n matrix to store to the file
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


void storeMatrix(const char *fileName, const int m, const int n, const double *matrix)
{

#ifdef PROFILE_STORE_MATRIX
	struct timeval start, end;
	gettimeofday(&start, 0);
#endif
#if DEBUG_LEVEL >= 2
	printf("Entering storeMatrix\n");
#endif
#ifdef ERROR_CHECKING
	errno = 0;		//no error occured so far
#endif

	int i, j;
	FILE *file = fopen(fileName, "w");

#ifdef ERROR_CHECKING
	if (errno) {
	  perror(fileName);
	  return;
	}
#endif

	fprintf(file, "%d\n", m);
	fprintf(file, "%d\n", n);

	for(i = 0; i < m; i++)
	{
		for(j = 0; j < n; j++)
		{
			fprintf(file, "%.15e ", matrix[i + j * m]);
		}
		fprintf(file, "\n");
	}
	fclose(file);

#if DEBUG_LEVEL >= 2
	printf("Exiting storeMatrix\n");
#endif
#ifdef PROFILE_STORE_MATRIX
	gettimeofday(&end, 0);
	outputTiming("Timing:", start, end);
#endif

}
//end of storeMatrix
//------------------
