#ifndef __FREQUENCY_ESTIMATOR_H
#define __FREQUENCY_ESTIMATOR_H

#include <vector>
#include <cstring>
#include <cmath>

#include "MathGenMin.h"
#include "Error.h"
#include "bcf_filtered_reader.h"
#include "Eigen/Dense"
#include "Eigen/Core"
#include "Eigen/SVD"

class frequency_estimator { //: public VectorFunc {
 public:

  Eigen::MatrixXd* pEigVecs;
  Eigen::BDCSVD<Eigen::MatrixXd>* pSVD;  
  
  int32_t nsamples;
  int32_t ndims;
  double tol;
  double maxLambda;
  
  bcf_hdr_t* hdr;
  bcf1_t* iv;
  float hwe0z;
  float hwe1z;
  float ibc0;
  float ibc1;

  int32_t* pls;
  int32_t n_pls;
  int8_t* ploidies;
  float* ifs;

  double pooled_af;
  bool isaf_computed;

  frequency_estimator(Eigen::MatrixXd* _pEigVecs, double _tol = 1e-10, double maxLambda = 1.0);
  frequency_estimator(Eigen::BDCSVD<Eigen::MatrixXd>* pSVD, double _tol = 1e-10, double maxLambda = 1.0);  
  ~frequency_estimator();

  bool set_hdr(bcf_hdr_t* _hdr);
  bool set_variant(bcf1_t* _iv, int8_t* ploidies, int32_t* _pl = NULL);
  double estimate_pooled_af_em(int32_t maxiter=10);
  
  double estimate_isaf_em(int32_t maxiter = 20);
  //double estimate_isaf_simplex();
  
  bool score_test_hwe(bool use_isaf = true);
  bool update_variant();
  //virtual double Evaluate(Vector& v);  

  /*

  double testHWE(Vector& v) {
    int32_t nsamples = pEst->pBfr->get_nsamples();
    double expGeno, isaf, isafQ, llk;
    
    double* smEigVecs;
    int32_t i, j;
    double p0, p1, p2, U, sumU, sumU2;
    BCFFilteredReader& bfr = *(pEst->pBfr);
    
    llk = 0; pEst->meanISAF = 0; pEst->maxISAF = 0; pEst->minISAF = 1;
    if ( pEst->ndims+1 != v.dim )
      error("[E:%s:%d %s] Dimensions do not match %d vs %d", __FILE__, __LINE__, __PRETTY_FUNCTION__, pEst->ndims, v.dim);
    
    sumU = sumU2 = 0;
    
    for(i=0; i < nsamples; ++i) {
      expGeno = v[0];
      
      smEigVecs = pEst->pEigVecs->at(i);
      for(j=0; j < pEst->ndims; ++j) {
	expGeno += (v[j+1] * smEigVecs[j]);
      }
      
      isaf = expGeno/2.0;
      //isaf = invLogit(expGeno/2.0);	
      if ( isaf < pEst->minMAF ) isaf = pEst->minMAF;
      else if ( 1.0-isaf < pEst->minMAF ) isaf = 1.0-pEst->minMAF;
      
      if ( pEst->maxISAF < isaf ) pEst->maxISAF = isaf;
      if ( pEst->minISAF > isaf ) pEst->minISAF = isaf;
      pEst->meanISAF += isaf;
      
      isafQ = 1.0-isaf;
      p0 = phredConv.toProb(bfr.get_likelihood_at(i*3));
      p1 = phredConv.toProb(bfr.get_likelihood_at(i*3+1));
      p2 = phredConv.toProb(bfr.get_likelihood_at(i*3+2));
      
      U = (p0-2*p1+p2)/(isafQ*isafQ*p0+2*isaf*isafQ*p1+isaf*isaf*p2+1e-100);
      sumU += U;
      sumU2 += (U*U);
    }
    pEst->meanISAF /= nsamples;      
    return sumU/sqrt(sumU2);
  }
  
  virtual double Evaluate(Vector& v) {
    int32_t nsamples = pEst->pBfr->get_nsamples();
    double expGeno, isaf, isafQ, llk;
    
    double* smEigVecs;
    int32_t i, j;
    BCFFilteredReader& bfr = *(pEst->pBfr);
    
    llk = 0;
    if ( pEst->ndims+1 != v.dim )
      error("[E:%s:%d %s] Dimensions do not match %d vs %d", __FILE__, __LINE__, __PRETTY_FUNCTION__, pEst->ndims, v.dim);
    
    for(i=0; i < nsamples; ++i) {
      expGeno = v[0];
      
      smEigVecs = pEst->pEigVecs->at(i);
      for(j=0; j < pEst->ndims; ++j) {
	expGeno += (v[j+1] * smEigVecs[j]);
      }
      
      isaf = expGeno/2.0;
      //isaf = invLogit(expGeno/2.0);
      if ( isaf < pEst->minMAF ) isaf = pEst->minMAF;
      else if ( 1.0-isaf < pEst->minMAF ) isaf = 1.0-pEst->minMAF;
      
      isafQ = 1.0-isaf;
      llk += log(isafQ * isafQ * phredConv.toProb(bfr.get_likelihood_at(i*3)) + 2.0 * isaf * isafQ * phredConv.toProb(bfr.get_likelihood_at(i*3+1)) + isaf * isaf * phredConv.toProb(bfr.get_likelihood_at(i*3+2)));
      }
      return 0-llk;
    }
  };
  */
};

#endif
