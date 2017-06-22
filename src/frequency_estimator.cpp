#include "frequency_estimator.h"

frequency_estimator::frequency_estimator(Eigen::MatrixXd* _pEigVecs, double _tol, double _maxLambda) {
  if ( _pEigVecs == NULL )
    error("[E:%s:%d %s] Invalid eigenvectors", __FILE__, __LINE__, __PRETTY_FUNCTION__ );

  pEigVecs = _pEigVecs;
  nsamples = (int32_t)pEigVecs->rows();
  ndims = (int32_t)pEigVecs->cols();

  pSVD = new Eigen::BDCSVD<Eigen::MatrixXd>(*_pEigVecs, Eigen::ComputeThinU | Eigen::ComputeThinV);

  // set dimensions;
  tol = _tol;
  maxLambda = _maxLambda;  
  
  hdr = NULL; iv = NULL;
  hwe0z = hwe1z = ibc0 = ibc1 = 0;

  pls = NULL; n_pls = 0; ploidies = NULL;
  ifs = new float[nsamples];
  pooled_af = -1;
  isaf_computed = false;
}

frequency_estimator::~frequency_estimator() {
  if ( ( pEigVecs ) && ( pSVD ) )
    delete pSVD;
  if (ifs != NULL )
    delete [] ifs;  
}

frequency_estimator::frequency_estimator(Eigen::BDCSVD<Eigen::MatrixXd>* _pSVD, double _tol, double _maxLambda) {
  // set dimensions;
  pEigVecs = NULL;

  pSVD = _pSVD;
  nsamples = (int32_t)pSVD->matrixU().rows();
  ndims = (int32_t)pSVD->matrixU().cols();

  tol = _tol;
  maxLambda = _maxLambda;    
  
  hdr = NULL; iv = NULL;
  hwe0z = hwe1z = ibc0 = ibc1 = 0;

  pls = NULL; n_pls = 0; ploidies = NULL;
  ifs = new float[nsamples];
  pooled_af = -1;
  isaf_computed = false;  
}

bool frequency_estimator::set_hdr(bcf_hdr_t* _hdr) {
  if ( hdr != _hdr ) {
    hdr = _hdr;
    char buffer[65535];
    sprintf(buffer,"##INFO=<ID=IS_HWE_SLP,Number=1,Type=Float,Description=\"Z-score of HWE test with individual-sepcific allele frequencyes\">\n");
    bcf_hdr_append(hdr, buffer);
    sprintf(buffer,"##INFO=<ID=IS_IBC,Number=1,Type=Float,Description=\"Inbreeding coefficient with individual-sepcific allele frequencies\">\n");
    bcf_hdr_append(hdr, buffer);
    sprintf(buffer,"##INFO=<ID=MAX_IF,Number=1,Type=Float,Description=\"Maximum Individual-specific allele frequency\">\n");
    bcf_hdr_append(hdr, buffer);    
    sprintf(buffer,"##INFO=<ID=MIN_IF,Number=1,Type=Float,Description=\"Minimum Individual-specific allele frequency\">\n");
    bcf_hdr_append(hdr, buffer);    
    sprintf(buffer,"##FORMAT=<ID=IF,Number=1,Type=Float,Description=\"Individual-specific allele frequencies\">\n");
    bcf_hdr_append(hdr, buffer);    
    bcf_hdr_sync(hdr);
    return true;
  }
  else return false;
}

bool frequency_estimator::set_variant(bcf1_t* _iv, int8_t* _ploidies, int32_t* _pl) {
  iv = _iv;
  ploidies = _ploidies;
  if ( iv->n_sample != nsamples )
    error("[E:%s:%d %s] nsamples %d != %d in the EigenVector", __FILE__, __LINE__, __PRETTY_FUNCTION__, iv->n_sample, nsamples);
    
  // parse PL fields
  bcf_unpack(iv, BCF_UN_ALL);
    
  if ( _pl != NULL ) { pls = _pl; n_pls = 3*nsamples; }
  else if ( bcf_get_format_int32(hdr, iv, "PL", &pls, &n_pls) < 0 ) {
    float gls[nsamples+1];
    if ( bcf_get_format_float(hdr, iv, "GL", &gls, &n_pls) < 0 ) {
      error("[E:%s:%d %s] Cannot parse PL or GL field", __FILE__, __LINE__, __PRETTY_FUNCTION__);
    }
    else {
      if ( pls == NULL ) pls = new int32_t[n_pls];      
      for(int32_t i=0; i < nsamples; ++i) {
	float maxgl = gls[3*i];
	if ( gls[3*i+1] > maxgl ) maxgl = gls[3*i+1];
	if ( gls[3*i+2] > maxgl ) maxgl = gls[3*i+2];
	pls[3*i] = (int32_t)floor(-10*(gls[3*i]-maxgl)+0.5);
	pls[3*i+1] = (int32_t)floor(-10*(gls[3*i+1]-maxgl)+0.5);
	pls[3*i+2] = (int32_t)floor(-10*(gls[3*i+2]-maxgl)+0.5);	
      }
    }
  }

  pooled_af = -1;
  isaf_computed = false;
  
  return true;
  //else return false;
}

double frequency_estimator::estimate_pooled_af_em(int32_t maxiter) {
  if ( pooled_af < 0 ) {
    double p = 0.5, q = 0.5;
    double p0, p1, p2, sump, apsum, an;
    for(int32_t i=0; i < maxiter; ++i) {
      apsum = 0;
      an = 0;
      for(int32_t j=0; j < nsamples; ++j) {
	if ( ploidies[j] == 2 ) {
	  p0 = q*q*phredConv.toProb(pls[j*3]);
	  p1 = 2*p*q*phredConv.toProb(pls[j*3+1]);
	  p2 = p*p*phredConv.toProb(pls[j*3+2]);
	  sump = p0+p1+p2;
	  apsum += (p1/sump + p2/sump*2.0);
	  an += 2;
	}
	else if ( ploidies[j] == 1 ) {
	  p0 = q*phredConv.toProb(pls[j*3]);
	  p2 = p*phredConv.toProb(pls[j*3+2]);
	  sump = p0+p2;
	  apsum += (p2/sump);
	  ++an;
	}
      }
      p = apsum/an;
      q = 1.0-p;
    }
    pooled_af = p;
    if ( pooled_af * nsamples < 0.5 ) { pooled_af = 0.5 / ( 1 + 2 *nsamples ); }
    else if ( (1-pooled_af)*nsamples < 0.5 ) { pooled_af = 1 - 0.5/(1 + 2*nsamples); }
  }
  return pooled_af;
}

bool frequency_estimator::score_test_hwe(bool use_isaf) {
  estimate_pooled_af_em();    
    
  double pp0 = (1.-pooled_af)*(1.-pooled_af);
  double pp1 = 2*pooled_af*(1-pooled_af);
  double pp2 = pooled_af*pooled_af;
  double sumU0 = 0, sqU0 = 0, sumU1 = 0, sqU1 = 0;
  double obsH0 = 0, obsH1 = 0, expH1 = 0;
  double l0, l1, l2, sum1, sum0, ph1, ph0, U0, U1;
  int32_t ndiploids = 0;

  // pretend that we have pp0, pp1, pp2 observations of each genotype (pseudocount)
  sumU0 = sumU1 = pp0*(1/pp0) + pp1 * (-2/pp1) + pp2*(1/pp2);
  sqU0 = sqU1 = pp0*1/pp0/pp0 + pp1 * 4/pp1/pp1 + pp2*(1/pp2/pp2);

  for(int32_t i=0; i < nsamples; ++i) {
    if ( ploidies[i] != 2 ) continue;
    ++ndiploids;
    l0 = phredConv.toProb(pls[i*3]);
    l1 = phredConv.toProb(pls[i*3+1]);
    l2 = phredConv.toProb(pls[i*3+2]);

    sum0 = l0*pp0 + l1*pp1 + l2*pp2 + 1e-100;
    ph0 = l1*pp1;
    obsH0 += (ph0/sum0);
    U0 = (l0-2*l1+l2)/sum0;
    sumU0 += U0;
    sqU0 += (U0*U0);
    
    if ( use_isaf ) {
      estimate_isaf_em();
      
      sum1 = l0*(1.-ifs[i])*(1.-ifs[i]) + 2*l1*(1.-ifs[i])*ifs[i] + l2*ifs[i]*ifs[i] + 1e-100;
      ph1 = 2*l1*(1.-ifs[i])*ifs[i];
      obsH1 += (ph1/sum1);
      expH1 += (2*(1.-ifs[i])*ifs[i]);
      U1 = (l0-2*l1+l2)/sum1;
      sumU1 += U1;
      sqU1 += (U1*U1);
    }
  }

  hwe0z = sumU0/sqrt(sqU0);
  ibc0 = 1.0 - (obsH0+1)/(pp1*ndiploids+1);

  if ( use_isaf ) {
    hwe1z = sumU1/sqrt(sqU1);
    ibc1 = 1.0 - (obsH1+1)/(expH1+1);
  }

  return true;
}


/*
double frequency_estimator::Evaluate(Vector& v) {
  double llk = 0;
  if ( ndims+1 != v.dim )
    error("[E:%s:%d %s] Dimensions do not match %d vs %d", __FILE__, __LINE__, __PRETTY_FUNCTION__, ndims, v.dim);

  double expGeno, isaf, isafQ;
    
  for(int32_t i=0; i < nsamples; ++i) {
    expGeno = v[0];
    for(int32_t j=0; j < ndims; ++j)
      expGeno += (pEigVecs->operator()(i,j) * v[j+1]);
      
    isaf = expGeno/2.0;
    if ( isaf*(nsamples+nsamples+1) < 0.5 ) isaf = 0.5/(nsamples+nsamples+1);
    if ( (1.0-isaf)*(nsamples+nsamples+1) < 0.5 ) isaf = 1.0-0.5/(nsamples+nsamples+1);
    isafQ = 1.0-isaf;

    ifs[i] = isaf;
    if ( ploidies[i] == 2 ) {
      llk += log(isafQ * isafQ    * phredConv.toProb(pls[i*3]) +
		 2 * isafQ * isaf * phredConv.toProb(pls[i*3+1]) +
		 isaf * isaf      * phredConv.toProb(pls[i*3+2]));
    }
    else if ( ploidies[i] == 1 ) {
      llk += log(isafQ * phredConv.toProb(pls[i*3]) +
		 isaf  * phredConv.toProb(pls[i*3+2]));
    }

  }
  return 0-llk;  
}
*/

double frequency_estimator::estimate_isaf_em(int32_t maxiter) {
  //Eigen::MatrixXd eV = Eigen::MatrixXd::Constant(nsamples,ndims+1,1.0);
  //eV.block(0,1,nsamples,ndims) = *pEigVecs;
  //Eigen::BDCSVD<Eigen::MatrixXd> svd(eV, Eigen::ComputeThinU | Eigen::ComputeThinV);
  if ( !isaf_computed ) {
    estimate_pooled_af_em();     
    Eigen::VectorXd y(nsamples);
    Eigen::VectorXd isaf = Eigen::VectorXd::Constant(nsamples, pooled_af);

    double maf = pooled_af > 0.5 ? 1-pooled_af : pooled_af;
    //Eigen::VectorXd lambda = Eigen::VectorXd::Constant(ndims, maxLambda * (1.-maf) / (maf * nsamples * 2));
    double lambda = maxLambda * (1.-maf) / (maf * nsamples * 2);
    double p0, p1, p2;

    for(int32_t i=0; i < maxiter; ++i) { // maxiter = 30
      for(int32_t j=0; j < nsamples; ++j) {
	if ( ploidies[j] == 2 ) {
	  p0 = ( 1.0 - isaf[j] ) * ( 1.0 - isaf[j] ) * phredConv.toProb(pls[3*j]);
	  p1 = 2 * isaf[j] * ( 1.0 - isaf[j] ) * phredConv.toProb(pls[3*j+1]);
	  p2 = isaf[j] * isaf[j] * phredConv.toProb(pls[3*j+2]);
	  y[j] = (p1+p2+p2+1e-100)/(p0+p1+p2+1e-100);
	}
	else if ( ploidies[j] == 1 ) {
	  p0 = ( 1.0 - isaf[j] ) * phredConv.toProb(pls[3*j]);
	  p2 = isaf[j] * phredConv.toProb(pls[3*j+2]);
	  y[j] = (p2+p2+1e-100)/(p0+p2+1e-100);
	}
	else {
	  y[j] = isaf[j] * 2;
	}
      }
      // U diag(d_i^2/(d_i^2+lambda)) U'y
      Eigen::VectorXd d2 = pSVD->singularValues();
      Eigen::MatrixXd UD2 = pSVD->matrixU(); //
      for(int32_t j=0; j < nsamples; ++j) {
	for(int32_t k=0; k < ndims; ++k) {
	  UD2(j,k) *= ( d2[k] / ( d2[k] + lambda ) );
	}
      }
      
      isaf = UD2 * ( pSVD->matrixU().transpose() * y ) / 2.0;
      for(int32_t j=0; j < nsamples; ++j) {
	if ( isaf[j]*(nsamples+nsamples+1) < 0.5 ) isaf[j] = 0.5/(nsamples+nsamples+1);
	else if ( (1.0-isaf[j])*(nsamples+nsamples+1) < 0.5 ) isaf[j] = 1.0-0.5/(nsamples+nsamples+1);
      }
    }

    for(int32_t j=0; j < nsamples; ++j) {
      ifs[j] = (float)isaf[j];
    }
    isaf_computed = true;
  }
  
  return 0;
}

/*
double frequency_estimator::estimate_isaf_simplex() {
  AmoebaMinimizer isafMinimizer;
  Vector startingPoint(ndims+1);
  double emaf = estimate_pooled_af_em();
  
  startingPoint.Zero();
  startingPoint[0] = emaf*2.0;

  isafMinimizer.func = this;

  isafMinimizer.Reset(ndims+1);
  isafMinimizer.point = startingPoint;
  isafMinimizer.Minimize(tol);

  score_test_hwe(true, emaf);

  //iter = isafMinimizer.cycleCount;
  return 0-isafMinimizer.fmin;
  //return ancMinimizer.fmin;
}
*/

 
bool frequency_estimator::update_variant() {
  float hweslp0 = (float)((hwe0z > 0 ? -1 : 1) * log10( erfc(fabs(hwe0z)/sqrt(2.0)) + 1e-100 ));
  float hweslp1 = (float)((hwe1z > 0 ? -1 : 1) * log10( erfc(fabs(hwe1z)/sqrt(2.0)) + 1e-100 ));
  float max_if = 0, min_if = 1;
  for(int32_t j=0; j < nsamples; ++j) {
    if ( ifs[j] > max_if ) max_if = ifs[j];
    if ( ifs[j] < min_if ) min_if = ifs[j];
  }
    
  bcf_update_info_float(hdr, iv, "HWE_SLP", &hweslp0, 1);
  bcf_update_info_float(hdr, iv, "IBC", &ibc0, 1);
  bcf_update_info_float(hdr, iv, "IS_HWE_SLP", &hweslp1, 1);
  bcf_update_info_float(hdr, iv, "IS_IBC", &ibc1, 1);
  bcf_update_info_float(hdr, iv, "MAX_IF", &max_if, 1);
  bcf_update_info_float(hdr, iv, "MIN_IF", &min_if, 1);    
  bcf_update_format_float(hdr, iv, "IF", ifs, nsamples);
  return true;
}

