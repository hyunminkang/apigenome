/** loadMatrix
 *
 * Purpose:
 * 		Loads a matrix from the specified file and allocates memory according to the matrix dimensions
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
 * fileName	in, 	The path to the file to read
 *
 * m		out, 	first dimension of the loaded matrix
 *
 * n		out, 	second dimension of the loaded matrix
 *
 * matrix	out, 	m*n matrix loaded from the file
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


void loadMatrix(const char *fileName, int *m, int *n, double **matrix)
{
#ifdef PROFILE_LOAD_MATRIX
	struct timeval start, end;
	gettimeofday(&start, 0);
#endif
#if DEBUG_LEVEL >= 2
	printf("Exiting loadMatrix\n");
#endif

#ifdef ERROR_CHECKING	
	errno = 0;		//no error occured so far
#endif

	int i, j, matches;
	FILE *file = fopen(fileName, "r");

#ifdef ERROR_CHECKING
	if (errno) {
	  perror(fileName);
	  return;
	}
#endif

	matches = fscanf(file, "%d", m);
	matches = fscanf(file, "%d", n);

	(*matrix) = (double*) malloc(sizeof(double) * (*m) * (*n));

#ifdef ERROR_CHECKING
	if (errno) {
	  perror("Allocating memory in loadMatrix failed");
	  return;
	}
#endif

	double tmp;
	for(i = 0; i < *m; i++)
	{
		for(j = 0; j < *n; j++)
		{
			matches = fscanf(file, "%lf", &tmp);
			(*matrix)[i + j * (*m)] = tmp;
		}
	}
	fclose(file);
#if DEBUG_LEVEL >= 2
	printf("Exiting loadMatrix\n");
#endif
#ifdef PROFILE_LOAD_MATRIX
	gettimeofday(&end, 0);
	outputTiming("Timing:", start, end);
#endif
}
//end of loadMatrix
//-----------------
