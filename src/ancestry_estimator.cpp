#include "ancestry_estimator.h"

ancestry_estimator::ancestry_estimator(std::vector<double*>* _pLoadings, std::vector< std::vector<double> >* _pProbs, int32_t _ndims, double _minMAF, double _tol) 
  : pLoadings(_pLoadings), pProbs(_pProbs), ndims(_ndims), minMAF(_minMAF), tol(_tol), llk0(0), llk1(0), iter(-1) {
  if ( _pLoadings == NULL ) error("[E:%s:%d %s] Invalid PC loadings", __FILE__, __LINE__, __PRETTY_FUNCTION__ );
  if ( _pProbs == NULL ) error("[E:%s:%d %s] Invalid genotype likelioods", __FILE__, __LINE__, __PRETTY_FUNCTION__ );
  if ( pLoadings->size() * 3 != pProbs->at(0).size() )
    error("Number of markers do not match between PC loadings and likelihoods : %u vs %u", pLoadings->size(), pProbs->at(0).size());
  //notice("[%s] Finishing",__PRETTY_FUNCTION__);
}

double ancestry_estimator::optimizeLLK(int32_t idx, double* optPC) {
  pc_llk_func lf(this, idx);
  AmoebaMinimizer ancMinimizer;
  Vector startingPoint(ndims);
  startingPoint.Zero();
  ancMinimizer.func = &lf;

  llk0 = 0-lf.Evaluate(startingPoint);
  
  ancMinimizer.Reset(ndims);
  ancMinimizer.point = startingPoint;
  ancMinimizer.Minimize(tol);

  /*
  std::string msg = "optimizeLLK";
  catprintf(msg, "(%d) LLK = %lg, PC = (", idx, ancMinimizer.fmin);
  for(int32_t i=0; i < ndims; ++i) {
    catprintf(msg, "%s%.5lf", i == 0 ? "" : ", ", ancMinimizer.point[i]);
  }
  catprintf(msg,"), cycleCount=%d, cycleMax=%d", ancMinimizer.cycleCount, ancMinimizer.cycleMax);
  notice("%s", msg.c_str());
  */

  if ( optPC != NULL ) {
    for(int32_t i=0; i < ndims; ++i)
      optPC[i] = ancMinimizer.point[i];
  }

  iter = ancMinimizer.cycleCount;
  return (llk1 = 0-ancMinimizer.fmin);
  //return ancMinimizer.fmin;
}
