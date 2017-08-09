#ifndef __CONTAM_ESTIMATOR_H
#define __CONTAM_ESTIMATOR_H

#include <vector>
#include <cstring>
#include <cmath>
#include <stdint.h>

#include "MathGenMin.h"
#include "PhredHelper.h"
#include "Error.h"
#include "var_dict.h"
#include "hts_utils.h"

struct vb_plp {
  double* ud;
  int32_t acs[3];
  std::vector< std::pair<uint8_t, uint8_t> >* albqs;
  //std::vector<uint8_t> als;
  //std::vector<uint8_t> bqs;
  vb_plp() : ud(NULL), albqs(NULL) { acs[0] = acs[1] = acs[2] = 0; }
  
  int32_t add_al_bq(uint8_t al, uint8_t bq) {
    if ( albqs == NULL )
      albqs = new std::vector< std::pair<uint8_t, uint8_t> >();
    albqs->push_back( std::pair<uint8_t,uint8_t>(al, bq) );
    ++(acs[al]);
    return (int32_t)albqs->size();
  }
  
  int32_t depth() const {
    return ( albqs == NULL ? 0 : (int32_t)albqs->size() );
  }

  double get_isaf(int32_t dim, const double* v, double minMAF) const {
    double expGeno = ud[0] * 2.0;
    //notice("ud = %lf", ud[0]);
    for(int32_t i=0; i < dim; ++i) {
      //notice("i=%d ud=%lf v=%lf", i, ud[i+1], v[i]);
      expGeno += (ud[i+1] * v[i]);
    }
    expGeno /= 2.0;
    if ( expGeno < minMAF ) return minMAF;
    else if ( expGeno + minMAF > 1 ) return 1-minMAF;
    else return expGeno;
  }

  // GLs
  double get_cont_gls(double alpha, double* gls, double* probRef) {
    int32_t i, j, k;
    int32_t d = albqs->size();
    uint8_t al, bq;
    double  err, mat, haf, sum, p[3];
    
    for(i=0; i < 9; ++i) gls[i] = 1.0;

    std::string als;
    std::string bqs;

    double logsum = 0;

    for(k=0; k < d; ++k) {
      al = albqs->at(k).first;
      bq = albqs->at(k).second;
      if ( al > 1 ) continue;

      als += (char)('0'+al);
      bqs += (char)(33+bq);

      err = phredConv.phred2Err[bq]/3.0;
      mat = phredConv.phred2Mat[bq];
      haf = 0.5 - err;

      // p[0] = Pr(Read | RR), p[1] = Pr(Read | RA), p[2] = Pr(Read | AA)
      if ( probRef == NULL ) {      
	p[0] = ( (al == 0) ? mat : err );
	p[1] = haf;
	p[2] = ( (al == 0) ? err : mat );
      }
      else {
	//error("foo");
	for(j=0; j < 3; ++j) {
	  p[j] = ( (al == 0) ? ( probRef[j] * mat + (1-probRef[j]) * err ) : ( probRef[j] * err + (1-probRef[j]) * mat ) );
	}
      }

      // Pr(R|G1,G2,alpha) = (1-alpha)Pr(R|G1) + alpha*Pr(R|G2)
      sum = 0;
      for(i=0; i < 3; ++i) { // intended genotype
	for(j=0; j < 3; ++j) {  
	  sum += (gls[i*3+j] *= ((1-alpha) * p[i] + alpha * p[j]));
	}
      }
      for(i=0; i < 9; ++i)
	gls[i] /= sum;
      
      logsum += log(sum);

      //notice("k = %d, sum = %lf, alpha = %lf", k, sum, alpha);
    }

    for(i=0; i < 3; ++i) { // intended genotype
      for(j=0; j < 3; ++j) {
	if ( gls[i*3+j] < 1e-25 ) gls[i*3+j] = 1e-25;
      }
    }

    return logsum;
      
    notice("%.3lg %.3lg %.3lg %.3lg %.3lg %.3lg %.3lg %.3lg %.3lg %s %s",gls[0],gls[2],gls[3],gls[4],gls[5],gls[6],gls[7],gls[8],alpha,als.c_str(),bqs.c_str());
    if ( rand() % 1000 == 0 ) abort();
  }
};

typedef struct vb_plp vb_plp_t;

class contam_estimator {
 public:
  std::vector< var_elem<vb_plp> > plps;
  int32_t numPC;
  double minMAF;
  std::vector<double> muPCs;
  std::vector<double> sdPCs;

  double* fixPC;
  double  fixAlpha;
  bool    refBias;
  //bool    withinAncestry;
  
  double tol;
  double llk0;
  double llk1;
  //double llk2;
  int32_t iter;
  double probRef[3];
  //int32_t iter2;

 contam_estimator(int32_t _numPC, double _minMAF) : numPC(_numPC), minMAF(_minMAF), muPCs(_numPC,0), sdPCs(_numPC, 1), fixPC(NULL), fixAlpha(-1), refBias(false), tol(1e-10), llk0(0), llk1(0), iter(0) {
    probRef[0] = 1.0;
    probRef[1] = 0.5;
    probRef[2] = 0.0;
  }

  int32_t remove_low_depth(int32_t minDP, int32_t* v_rmv = NULL,
			   int32_t* v_kep = NULL, int32_t* r_rmv = NULL,
			   int32_t* r_kep = NULL) {
    int32_t nv_rmv = 0, nv_kep = 0, nr_rmv = 0, nr_kep = 0;
    std::vector< var_elem<vb_plp> > new_plps;
    for(int32_t i=0; i < (int32_t)plps.size(); ++i) {
      if ( plps[i].val.depth() >= minDP ) {
	nr_kep += plps[i].val.depth();
	++nv_kep;
	new_plps.push_back(plps[i]);
      }
      else {
	// remove allocated memory
	nr_rmv += plps[i].val.depth();
	++nv_rmv;
	delete [] plps[i].val.ud;
	delete plps[i].val.albqs;
      }
    }
    plps.swap(new_plps);
    
    if ( v_rmv ) *v_rmv = nv_rmv;
    if ( v_kep ) *v_kep = nv_kep;
    if ( r_rmv ) *r_rmv = nr_rmv;
    if ( r_kep ) *r_kep = nr_kep;    
    
    return nv_rmv;
  }

  double optimizeLLK(double* optAlpha, double* optPC1, double* optPC2, bool withinAncestry, bool refBiasFlag = false);

  double* get_relative_depths(double* optPC);

  void writePileup(htsFile* wf, double alpha, double* pc1, double* pc2);
  
  class cont_llk_func : public VectorFunc {
  public:
    contam_estimator* pEst;
    bool withinAncestry;
    
    cont_llk_func(contam_estimator& p) : pEst(&p), withinAncestry(false) {}
    
    virtual double Evaluate(Vector& v) {
      double  gls[9], p1[3], p2[3], p, q, lk;
      int32_t i, j, k;      
      int32_t nvars = (int32_t)pEst->plps.size();
      double  alpha = ( pEst->fixAlpha < 0 ) ? 1/(1+exp(0-v[0])) : pEst->fixAlpha;
      //alpha = 1/(1+exp(0-alpha));
      if ( alpha > 0.5 ) alpha = 1.0 - alpha;
      double* pc1   = ( pEst->fixPC == NULL ) ? (pEst->fixAlpha < 0 ? v.data+1 : v.data ) : pEst->fixPC;
      double* pc2   = ( withinAncestry || pEst->fixPC ) ? pc1 :
	( pEst->fixAlpha < 0 ? v.data + pEst->numPC + 1 : v.data + pEst->numPC );
      
      double llk = 0;
      int32_t nb0 = 0, nb1 = 0, nb2 = 0;
      //notice("nvars = %d", nvars);
      for(i=0; i < nvars; ++i) {
	vb_plp_t& plp = pEst->plps[i].val;
	llk += plp.get_cont_gls(alpha, gls, pEst->refBias ? pEst->probRef : NULL);

	//notice("gls = (%lf %lf %lf %lf %lf %lf %lf %lf %lf)", gls[0], gls[1], gls[2], gls[3], gls[4], gls[5], gls[6], gls[7], gls[8]);
	
	p = plp.get_isaf(pEst->numPC, pc1, pEst->minMAF);
	q = 1-p;

	if ( ( plp.ud[0] <= pEst->minMAF ) || ( 1-plp.ud[0] <= pEst->minMAF ) ) ++nb0;

	if ( ( p <= pEst->minMAF ) || ( q <= pEst->minMAF ) ) {	++nb1; }

	p1[0] = q*q; p1[1] = 2*p*q; p1[2] = p*p;
	if ( pc1 == pc2 ) {
	  p2[0] = p1[0]; p2[1] = p1[1]; p2[2] = p1[2];
	}
	else {
	  p = plp.get_isaf(pEst->numPC, pc2, pEst->minMAF);
	  q = 1-p;

	  if ( ( p <= pEst->minMAF ) || ( q <= pEst->minMAF ) ) ++nb2;
	  
	  p2[0] = q*q; p2[1] = 2*p*q; p2[2] = p*p;	  
	}

	lk = 0;
	for(j=0; j < 3; ++j) {
	  for(k=0; k < 3; ++k) {
	    lk += (p1[j] * p2[k] * gls[j*3+k]);
	  }
	}
	//error("i = %d, lk = %lf, p = %lf, q = %lf", i, lk, p, q);
	llk += log(lk);
      }

      //notice("pc1 = (%.5lg, %.5lg), pc2 = (%.5lg, %.5lg), nb = (%d, %d, %d) / %d", pc1[0], pc1[1], pc2[0], pc2[1], nb0, nb1, nb2, nvars);
      //if ( nb1 + nb2 > 0 )
      notice("Evaluate: alpha = (%.5lf), pc1 = (%.5lg %.5lg), pc2 = (%.5lg %.5lg), llk = %.5lg, dim = %d, %d, %d, %d", alpha, pc1[0], pc1[1], pc2[0], pc2[1], llk, v.dim, nb0, nb1, nb2);	
	//if ( ( nb1 + nb2 ) > 2*nb0 + 0.05*nvars ) return 1e300;
      //if ( ( nb1 + nb2 > nvar ) || ( nb1 > 
	
      //notice("Evaluate: alpha = (%.5lf), pc1 = (%.5lg %.5lg), pc2 = (%.5lg %.5lg), llk = %.5lg, dim = %d", alpha, pc1[0], pc1[1], pc2[0], pc2[1], llk, v.dim);
      return 0-llk;
      //return llk;
    }
  };
};

#endif
