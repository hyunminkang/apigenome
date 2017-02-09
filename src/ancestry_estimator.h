#ifndef __ANCESTRY_ESTIMATOR_H
#define __ANCESTRY_ESTIMATOR_H

#include <vector>
#include <cstring>
#include <cmath>

#include "MathGenMin.h"
#include "Error.h"

class ancestry_estimator {
 public:
  std::vector<double*>* pLoadings;
  std::vector< std::vector<double> >* pProbs;
  int32_t ndims;
  double minMAF;
  double tol;
  double llk0;
  double llk1;
  int32_t iter;

  ancestry_estimator(std::vector<double*>* _pLoadings, std::vector< std::vector<double> >* _pProbs, int32_t _ndims, double _minMAF, double tol = 1e-10);
  double optimizeLLK(int32_t idx, double* optPC);
  
  class pc_llk_func : public VectorFunc {
  public:
    ancestry_estimator* pEst;
    int32_t smIdx;
    pc_llk_func(ancestry_estimator* p, int32_t i) : pEst(p), smIdx(i) {
      if ( pEst == NULL )
	error("[E:%s:%d %s] Cannot specify NULL ancestry estimator", __FILE__, __LINE__, __PRETTY_FUNCTION__);
    }
    
    virtual double Evaluate(Vector& v) {
      int32_t nvars = (int32_t)pEst->pLoadings->size();
      std::vector<double>& vProbs = pEst->pProbs->at(smIdx);
      double expGeno, isaf, isafQ, llk;
      double* varLoadings;
      int32_t i, j;

      llk = 0;
      if ( pEst->ndims != v.dim )
	error("[E:%s:%d %s] Dimensions do not match %d vs %d", __FILE__, __LINE__, __PRETTY_FUNCTION__, pEst->ndims, v.dim);

      for(i=0; i < nvars; ++i) {
	varLoadings = pEst->pLoadings->at(i);
	expGeno = 2*varLoadings[0];
	for(j=0; j < pEst->ndims; ++j) 
	  expGeno += ( varLoadings[j+1] * v[j] );

	isaf = expGeno/2.0;
	if ( isaf < pEst->minMAF ) isaf = pEst->minMAF;
	else if ( 1.0-isaf < pEst->minMAF ) isaf = 1.0-pEst->minMAF;
	
	isafQ = 1.0-isaf;

	llk += log(isafQ * isafQ * vProbs[i*3] + 2.0*isaf*isafQ * vProbs[i*3+1] + isaf*isaf*vProbs[i*3+2]);
      }

      // temporary printing
      /*
      std::string msg = __PRETTY_FUNCTION__;
      catprintf(msg, "(");
      for(int32_t i=0; i < pEst->ndims; ++i) {
	catprintf(msg, "%s%.5lf", i == 0 ? "" : ", ", v[i]);
      }
      catprintf(msg,") = %lf", llk);
      notice("%s", msg.c_str());*/
      /////////////////////////////////
      
      return 0-llk;
    }
  };
};

#endif
