/** outputTiming - Outputs profile timings
 *
 * Purpose:
 *		Calculate time difference between two stopped times and print the resulting timing
 *
 * Description:
 *		The time is printed in microseconds
 *
 * Arguments:
 *
 * name		in, 	Text preceeding timing results
 *
 * start	in, 	start of the timing interval
 *
 *end		in, 	end of the timing interval
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>


void outputTiming(const char *name, const struct timeval start, struct timeval end)
{
	struct timeval diff;
	diff.tv_sec = end.tv_sec - start.tv_sec;

	diff.tv_usec = end.tv_usec - start.tv_usec;

	printf("%d \t", (int)(((diff.tv_sec*1000000)+diff.tv_usec)));
}
//end of outputTiming
//-------------------
