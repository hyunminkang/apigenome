/**
 * randnumber - returns a random number between min and max
 *
 * Purpose:
 *		Creating random numbers for initialisation of factor matrices W and H
 *
 * Description:
 *		This routine uses the random number generator rand() from the c standard library
 * 
 * Arguments:
 * 
 * min			in, 	lower interval bound for random numbers
 *
 * max			in, 	upper interval bound for random numbers
 *
 * return value		out, 	generated random number
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>


double randnumber(const int min, const int max) 
{
  static int initialised = 0;			
  if (!initialised) {					
    srand(time(NULL));
    initialised = 1;
  }
  return (min + ( (max - min + 1.) * rand() / (double)(RAND_MAX) ));
}
//end of randnumber
//----------------
