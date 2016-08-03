//example_withoutdriver.c
//version 1.03
//#########################################
//#				  	  #
//#  Calling the computational routines	  #
//#				          #
//#########################################


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
#include <stdlib.h>
#include "loadmatrix.h"
#include "nmf_als.h"
#include "nmf_mu.h"


/** main function
 *
 */
int main(int argc, char** argv) {

  
// @ 1. Step: creating necessary data structures
// #############################################
  
  
  //matrix dimensions
  int m, n, k;

  //pointer to matrices
  double *A, *W0, *H0;

  //probably a good idea to check for conforming matrix dimensions
  
  loadMatrix("emaildata.matrix", &m, &n, &A);				// loads matrix "a.dat" into "A" and rows(A) = m, cols(A) = n
  loadMatrix("w0_rand_64.matrix", &m, &k, &W0);			// loads matrix "w0.dat" into "W0" and rows(W0) = m, cols(W0) = k
  loadMatrix("h0_rand_64.matrix", &k, &n, &H0);			// loads matrix "h0.dat" into "H0" and rows(H0) = k, cols(H0) = n

  // parameter for nmf computational routines
  int maxiter = 100;
  double tolerance = 2E-02;
  
// @ 2. Step: calling computational routines
//##########################################

  // @ 2.1: calculating NMF using MU and ALS algorithms
  //##########################################
  
  // calling parameters:		exemplary values		comment
  //##########################################################################################################################################
  
  // matrix A			"a.matrix"			Matrix to be factorized, A = W*H
  // matrix W0			W0				initialization factor matrix W0
  // matrix H0			H0				inititalization factor matrix H0
  // m				99				rows(A)
  // n				100				cols(A)
  // k				25				approximation factor
  // maxiter			100				maximum number of iterations
  // tolerance			2E-02				tolerance value for maxchange criterion
  // tolerance			2E-02				tolerance value for norm criterion



  // load and factorize matrix "a.dat" using random initialization and store resulting factor matrices
  //----------------------------------------------------------------------------------------------------
  
  // factor matrices are returned in W0 and H0 respectively
  // "maxiter" returns the number of actually run iterations
  // return value of every NMF computational routine is the frobenius norm of the result, scaled by dividing by (number of elements in A)

  double mu_norm = nmf_mu(A, &W0, &H0, m, n, k, &maxiter, tolerance, tolerance);
  printf("mu_norm = %e \t after iteration # %d\n", mu_norm, maxiter);
  
  
  
  // maxiter has to be reset
  maxiter = 100;
  
  
    // load and factorize matrix "a.dat" using random initialization and store resulting factor matrices
  //----------------------------------------------------------------------------------------------------
  
  // factor matrices are returned in W0 and H0 respectively
  // "maxiter" returns the number of actually run iterations
  // return value of every NMF computational routine is the frobenius norm of the result, scaled by dividing by (number of elements in A)
  
  double als_norm = nmf_als(&A, &W0, &H0, m, n, k, &maxiter, tolerance, tolerance);
  printf("als_norm = %e \t after iteration # %d\n", als_norm, maxiter);

  
  
  free(A);
  free(W0);
  free(H0);

  return 0;
}
