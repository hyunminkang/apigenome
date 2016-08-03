// example.c
// version 1.03
//#################################
//#				  #
//# Calling the driver routine	  #
//#				  #
//#################################


// Compiling
//---------

// gcc -Wall example.c -c -I../include/

// Linking
//--------

// Linking requires the library libnnmf.a (or libnnmf.so)
// This library can be created using the provided Makefile and target 'lib' (or target 'shared' for libnnmf.so)

// Using generic blas/lapack routines:		-lnnmf -larpack -llapack -lblas(depending on system maybe -lgfortran)
// Using atlas blas/lapack routines:		-lnnmf -larpack -llapack -lf77blas -lcblas -latlas (depending on system maybe -lgfortran)
// Using goto blas and generic lapack		-lnnmf -larpack -llapack -lgoto(depending on system maybe -lgfortran)

#include <stdio.h>
#include <string.h>
#include "nmfdriver.h"
#include "loadmatrix.h"


/** main function
 *
 */
int main(int argc, char** argv) {
// @ 1. call - default options
//############################
  
  // calling parameters:		exemplary values		comment
  //##########################################################################################################################################
  
  // filename matrix a			"a.matrix"			Matrix to be factorized			
  // factor k				25				approximation rank, has to be smaller than original matrix dimensions
  // max. iterations			100				
  // filename initial matrix w0		NULL				set to NULL to generate a matrix in the first repetition as well
  // filename initial matrix h0		NULL				set to NULL to generate a matrix in the first repetition as well
  // algorithm				"als" 	 			set to one of the implemented algorithms
  // options to be used			NULL				set to NULL to use default options
  
  
  // load and factorize matrix "a.dat" using random initialization and store resulting factor matrices
  if ( ( argc > 8 ) || ( argc < 3 ) ) {
	printf("Usage: %s [input_matrix.file] [num.factor] [output_w_file = final_w.matrix] [output_h_file = final_h.matrix] [algorithm=alspg,mu,als,neals,pg,fast,bayes] [init=rand,nndsvd,infogain,gainratio] [max.iteration=100] [initial.w0=NULL] [initial.h0=NULL]\n",argv[0]);
	return -1;
  }

  alg_t alg;
  if ( argc > 5 ) {
    if ( strcmp(argv[5],"mu") == 0 ) alg = mu;
    else if ( strcmp(argv[5],"als") == 0 ) alg = als;
    else if ( strcmp(argv[5],"neals") == 0 ) alg = neals;
    else if ( strcmp(argv[5],"alspg") == 0 ) alg = alspg;
    else if ( strcmp(argv[5],"pg") == 0 ) alg = pg;
    else if ( strcmp(argv[5],"fast") == 0 ) alg = fast;
    else if ( strcmp(argv[5],"bayes") == 0 ) alg= bayes;
    else {
      fprintf(stderr,"Cannot recognize the algorithm %s\n",argv[5]); 
      return -1;
    }
    printf("Algorithm = %s\n",argv[5]);
  }
  else {
    alg = alspg;
    printf("Algorithm = alspg (default)\n");
  }

  options_t opts;
  opts.rep = 1;
  if ( argc > 6 ) {
    if ( strcmp(argv[6],"rand") == 0 ) opts.init = init_rand;
    else if ( strcmp(argv[6],"nndsvd") == 0 ) {
      opts.init = init_nndsvd;
      opts.nndsvd_maxiter = -1;			//if set to -1 - default value will be set in generateMatrix
      opts.nndsvd_blocksize = 1;
      opts.nndsvd_tol = 2E-08;
      opts.nndsvd_ncv = -1;		
    }
    else if ( strcmp(argv[6],"infogain") == 0 ) opts.init = init_infogain;
    else if ( strcmp(argv[6],"gainratio") == 0 ) opts.init = init_gainratio;
    else {
      fprintf(stderr,"Cannot recognize the initialization method %s\n",argv[6]); 
      return -1;
    }
    printf("Initialization = init_%s\n",argv[6]);
  }
  else {
    opts.init = init_nndsvd;
    printf("Initialization = init_nndsvd (default)\n");
  }
  opts.min_init = 0;
  opts.max_init = 2;
  opts.w_out = (argc > 3) ? argv[3] : "final_w.matrix";
  opts.h_out = (argc > 4) ? argv[4] : "final_h.matrix";
  opts.TolX = 2.0E-06;
  opts.TolFun = 2.0E-06;

  nmfDriver(argv[1], atoi(argv[2]), (argc > 7) ? atoi(argv[7]) : 1000, (argc > 8) ? argv[8] : NULL, (argc > 9) ? argv[9] : NULL, alspg, &opts);
  /*  
// @ 2. call - explicitly set options: nndsvd initialization
//####################################
  
  
  // options is of type options_t* with following elements
  //######################################################
  
  // int rep;				1				Number of repetitions with new starting matrices
  // init_t init;			nndsvd				Method to use for initialising the starting matrices
  // int min_init;			0				minimal value for random initialisation
  // int max_init;			1				maximal value for random initialisaton
  // char* w_out;			"final_w.matrix"		Filename to write final matrix w to
  // char* h_out;			"final_h.matrix"		Filename to write final matrix h to
  // double TolX;			2.0E-02				tolerance value for convergence check of maxChange
  // double TolFun;			2.0E-02				tolerance value for convergence check of dnorm
  // int nndsvd_maxiter			-1				maxiter in nndsvd initialization (-1 requests setting of default value)
  // int nndsvd_blocksize		1				blocksize in nndsvd initialization; !! only works for 1 so far
  // double nndsvd_tol			2E-08				tolerance value for nndsvd initialization
  // int nndsvd_ncv			-1				length of arnoldi iteration (-1 requests setting of default value)
  
  // Creating option structure
  opts.rep = 1;
  opts.init = init_rand;
  opts.min_init = 0;
  opts.max_init = 2;
  opts.w_out = (argc >= 4) ? argv[3] : "final_w.matrix";
  opts.h_out = (argc >= 5) ? argv[4] : "final_h.matrix";
  opts.TolX = 2.0E-02;
  opts.TolFun = 2.0E-02;

  // load and factorize matrix "a.dat" using nndsvd initialization and store resulting factor matrices
  printf("ALSPG Algorithm, Default-Random Initialization:\n");
  nmfDriver(argv[1], atoi(argv[2]), (argc >= 6) ? atoi(argv[5]) : 100, (argc >=7) ? argv[6] : NULL, (argc >=8) ? argv[7] : NULL, mu, &opts);
  //nmfDriver("emaildata.matrix", 25, 100, NULL, NULL, mu, &opts);
  */
/*
  
  // @ 3. call - explicitly set options: gainratio initialization
//####################################
  
  
  // options is of type options_t* with following elements
  //######################################################
  
  // int rep;				1				Number of repetitions with new starting matrices
  // init_t init;			nndsvd				Method to use for initialising the starting matrices
  // int min_init;			0				minimal value for random initialisation
  // int max_init;			1				maximal value for random initialisaton
  // char* w_out;			"final_w.matrix"		Filename to write final matrix w to
  // char* h_out;			"final_h.matrix"		Filename to write final matrix h to
  // double TolX;			2.0E-02				tolerance value for convergence check of maxChange
  // double TolFun;			2.0E-02				tolerance value for convergence check of dnorm
  // int nndsvd_maxiter			-1				maxiter in nndsvd initialization (-1 requests setting of default value)
  // int nndsvd_blocksize		1				blocksize in nndsvd initialization; !! only works for 1 so far
  // double nndsvd_tol			2E-08				tolerance value for nndsvd initialization
  // int nndsvd_ncv			-1				length of arnoldi iteration (-1 requests setting of default value)
  
  // Creating option structure
  opts.rep = 1;
  opts.init = init_gainratio;
  opts.min_init = 0;
  opts.max_init = 1;
  opts.w_out = (argc >= 4) ? argv[3] : "final_w.matrix";
  opts.h_out = (argc >= 5) ? argv[4] : "final_h.matrix";
  opts.TolX = 2.0E-02;
  opts.TolFun = 2.0E-02;
  
  //load necessary class Column vector
  int m, n;
  loadMatrix("emailclasscolumn.matrix", &m, &n, &(opts.classColumn) );	// loads matrix "emailclasscolumn.matrix" into "opts.classColumn"
  opts.num_dist_classes = 3;						//3 mail classes -- alternative: count distinct classes in class column vector
  opts.parallelization = par_openmp;					//OpenMP parallelization || Default is pthreads


  // load and factorize matrix "a.dat" using nndsvd initialization and store resulting factor matrices
  printf("NEALS Algorithm, GainRatio Initialization using OpenMP parallelization:\n");  
  //nmfDriver("emaildata.matrix", 25, 100, NULL, NULL, neals, &opts);
  nmfDriver(argv[1], atoi(argv[2]), (argc >= 6) ? atoi(argv[5]) : 100, (argc >=7) ? argv[6] : NULL, (argc >=8) ? argv[7] : NULL, neals, &opts);
*/
  
  
  //free class Column vector.




  return 0;
}
