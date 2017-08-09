#include "contam_estimator.h"

void contam_estimator::writePileup(htsFile* wf, double alpha, double* pc1, double* pc2) {
  int32_t nvars = (int32_t)plps.size();
  hprintf(wf,"#CHROM\tPOS\tREF\tALT\tAF\tISAF1\tISAF2\tDEPTH\tALLELES\tBQ\n");
  for(int32_t i=0; i < nvars; ++i) {
    vb_plp& plp = plps[i].val;
    double af = plp.ud[0];
    double isaf1 = plp.get_isaf(numPC, pc1, minMAF);
    double isaf2 = plp.get_isaf(numPC, pc2, minMAF);
    int32_t d = plp.depth();
    std::string als, bqs;
    for(int32_t j=0; j < d; ++j) {
      als += (char)('0' + plp.albqs->at(j).first);
      bqs += (char)(33  + plp.albqs->at(j).second);
    }
    hprintf(wf,"%s\t%ld\t%s\t%s\t%.4lg\t%.4lg\t%.4lg\t%d\t%s\t%s\n",
	    plps[i].chr.c_str(), plps[i].var.pos1, plps[i].var.ref.c_str(),
	    plps[i].var.alts[0].c_str(), af, isaf1, isaf2, d,
	    als.c_str(), bqs.c_str());
  }
}

double* contam_estimator::get_relative_depths(double* optPC) {
  // Relative DP = [OBS_DP_AT_G]/[OBS_DP_AT_R]
  // With no sequencing error, we should expect
  // [E #R] = \sum_i d_i P_i(RR) + 0.5 * \sum_i d_i P_i(RA)
  // [E #A] = \sum_i d_i P_i(AA) + 0.5 * \sum_i d_i P_i(RA)
  // [O #R] 
  
  double ob[6] = {0, 0, 0, 0, 0, 0};
  double ex[6] = {0, 0, 0, 0, 0, 0};
  double gf[3], p, q;
  int32_t i, j;
  for(i=0; i < (int32_t)plps.size(); ++i) {
    vb_plp_t& plp = plps[i].val;
    p = plp.get_isaf(numPC, optPC, minMAF);
    q = 1-p;
    gf[0] = q*q; gf[1] = 2*p*q; gf[2] = p*p;    

    for(j=0; j < 3; ++j) {
      ob[j*2]   += (gf[j] * plp.acs[0]);
      ob[j*2+1] += (gf[j] * plp.acs[1]);
    }
    ex[0] += (gf[0] * (plp.acs[0]+plp.acs[1]));
    ex[2] += (gf[1] * (plp.acs[0]+plp.acs[1]) * 0.5);
    ex[3] += (gf[1] * (plp.acs[0]+plp.acs[1]) * 0.5);
    ex[5] += (gf[2] * (plp.acs[0]+plp.acs[1]));    
  }
  double* rdps = new double[15];
  for(i=0; i < 6; ++i) {
    rdps[3+i] = ob[i];
    rdps[9+i] = ex[i];    
  }
  for(i=0; i < 3; ++i) {
    rdps[i] = ( ob[i*2]/(ob[i*2]+ob[i*2+1]) ) / ( ex[i*2]/(ex[i*2]+ex[i*2+1]) );
  }
  return rdps;
}

double contam_estimator::optimizeLLK(double* optAlpha, double* optPC1, double* optPC2, bool withinAncestry, bool refBiasFlag) {
  refBias = refBiasFlag;
  
  int32_t ndims = ( fixAlpha < 0 ? 1 : 0 ) + ( ( fixPC == NULL ) ? ( withinAncestry ? numPC : 2*numPC ) : 0 ) + ( refBias ? 3 : 0 );

  if ( ndims == 0 ) {
    *optAlpha = fixAlpha;
    for(int32_t i=0; i < numPC; ++i)
      optPC1[i] = optPC2[i] = fixPC[i];
  }

  // calculate 

  Vector startingPoint(ndims);
  cont_llk_func lf(*this);
  
  if ( withinAncestry ) {
    startingPoint.Zero();
    if ( fixAlpha < 0 )
      startingPoint[0] = log(0.001/0.999);
    lf.withinAncestry = true;
  }
  else {
    int32_t o = 0;
    if ( fixAlpha < 0 ) {
      startingPoint[0] = log(*optAlpha/(1 - *optAlpha));
      o = 1;
    }
    for(int32_t i=0; i < numPC; ++i) {
      startingPoint[i+o] = (optPC1[i] - muPCs[i])/sdPCs[i];
      startingPoint[i+numPC+o] = (optPC2[i] - muPCs[i])/sdPCs[i];
    }
    lf.withinAncestry = false;
  }

  double epsilon = 1e-6;  
  if ( refBias ) {
    startingPoint[ndims-3] = log( (probRef[0]-epsilon) / (1-probRef[0]+epsilon) );
    startingPoint[ndims-2] = 0;
    startingPoint[ndims-1] = log( (probRef[2]+epsilon) / (1-probRef[2]-epsilon) );
  }

  AmoebaMinimizer conMinimizer;  
  conMinimizer.func = &lf;    

  llk0 = 0-lf.Evaluate(startingPoint);
  
  conMinimizer.Reset(ndims);
  conMinimizer.point = startingPoint;
  conMinimizer.Minimize(tol);

  *optAlpha = ( fixAlpha < 0 ) ? (1.0/(1+exp(0-conMinimizer.point[0]))) : fixAlpha;
  
  int32_t o1 = 0, o2 = 0;
  if ( *optAlpha > 0.5 ) {
    *optAlpha = 1 - *optAlpha;
    o1 = withinAncestry ? 0 : numPC;
    o2 = 0;
  }
  else {
    o1 = 0;        
    o2 = withinAncestry ? 0 : numPC;
  }
  
  if ( fixPC == NULL ) {
    for(int32_t i=0; i < numPC; ++i) {
      optPC1[i] = conMinimizer.point[i + (fixAlpha < 0 ? 1 : 0) + o1] * sdPCs[i] + muPCs[i];
      optPC2[i] = conMinimizer.point[i + (fixAlpha < 0 ? 1 : 0) + o2] * sdPCs[i] + muPCs[i];	
    }
  }

  if ( refBias ) {
    probRef[0] = 1.0 / (1.0 + exp(0-conMinimizer.point[ndims-3])) + epsilon;
    probRef[1] = 1.0 / (1.0 + exp(0-conMinimizer.point[ndims-2]));
    probRef[2] = 1.0 / (1.0 + exp(0-conMinimizer.point[ndims-1])) - epsilon;
    if ( probRef[0] > 1 ) probRef[0] = 1.0;
    if ( probRef[2] < 0 ) probRef[2] = 0.0;
  }
  
  iter = conMinimizer.cycleCount;
  llk1 = 0-conMinimizer.fmin;

  notice("llk0 = %.3lf, llk1 = %.3lf, iter = %d, alpha = %.3lg, refBias = (%.3lg, %.3lg, %.3lg)", llk0, llk1, iter, *optAlpha, probRef[0], probRef[1], probRef[2] );
  for(int32_t i=0; i < numPC; ++i) {
    notice("optPC1[%d] = %.5lg, optPC2[%d] = %.5lg", i, optPC1[i], i, optPC2[i]);
  }
  return llk1;
}
