/**
 * set_default_opts - Sets the default options in an options_t structure
 *
 * Purpose:
 *		If a null-pointer is passed to nmfDriver for the options-structure default values are used
 *		wich this routine implements
 * 
 * Description:
 *		The options structure contains following elements:
 *		rep		number of repetitions to run with differently initialised matrices	1
 *		init		method to use for initialising matrices					ran
 *		min_init	minimal value for random numbers in initialised matrices		0
 *		max_init	maximal value for random numbers in initialised matrices		1
 *		w_out		filename to store final matrix w in					"final_w.matrix"
 *		h_out		filename to store final matrix h in					"final_h.matrix"
 *		TolX		tolerance limit for maxChange						1E-04
 *		TolFun		tolerance for root mean square residual					1E-04
 *		nndsvd_maxiter	maximum iterations for SVD in ddsvd initialisation 			-1 -> default value set in generateMatrix
 *		nnsvd_blocksize	blocksize for SVD in ddsvd initialisation				64
 *		nndsvd_tol	tolerance for SVD in ddsvd initialisation				2E-16
 *		nndsvd_ncv	largest number of basis vectors in the Arnoldi Process			-1 -> ncv = 2 * nev (will be set in generateMatrix)
 *
 * Arguments:
 *
 * opts		in/out, 	option structure to store the default options in
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>

#include "common.h"


void set_default_opts_singleprec(options_t * opts)
{
  opts->rep = 1;
  opts->init = init_rand;
  opts->min_init = 0;
  opts->max_init = 1;
  opts->w_out = "final_w.matrix";
  opts->h_out = "final_h.matrix";
  opts->TolX = 1.0E-04f;
  opts->TolFun = 1.0E-04f;
  opts->nndsvd_maxiter = -1;			//if set to -1 - default value will be set in generateMatrix
  opts->nndsvd_blocksize = 64;
  opts->nndsvd_tol = 2E-16f;
  opts->nndsvd_ncv = -1;			//if set to -1 - default value will be set in generateMatrix
  opts->num_dist_classes = -1;
  opts->classColumn = NULL;
  opts->parallelization = par_pthreads;			//Method to use for parallelization: par_pthreads, par_openmp  
}
//end of set_default_opts
//-----------------------

