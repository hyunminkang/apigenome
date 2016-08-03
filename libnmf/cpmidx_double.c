/** cpmidx_double - compares two parameters of type idx_double, used for qsort sorting in descending order
 *
 * Purpose:
 *		Defines the comparison function needed by  c standard libraries "qsort()".
 *
 * Arguments:
 *
 * p1			in, 	first parameter to compare
 *
 * p2			in, 	second parameter to compare
 *
 * return value		out, 	returns -1 if p1 < p2, 1 if p1 > p2 and 0 if p1 = p2
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "common.h"


int cmpidx_double(const void *p1, const void *p2)
{
  idx_double *id1 = (idx_double*) p1;
  idx_double *id2 = (idx_double*) p2;

  if (id1->val < id2->val)
    return 1;
  else
    if (id1->val > id2->val)
      return -1;
    else return 0;
}
//end of cpmidx_double
//--------------------

