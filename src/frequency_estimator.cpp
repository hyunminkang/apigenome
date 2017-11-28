#include "frequency_estimator.h"

frequency_estimator::frequency_estimator(Eigen::MatrixXd* _pEigVecs, double _tol, double _maxLambda) :
  skipIf(false), skipInfo(false), siteOnly(false), gtError(0.005) 
{
  if ( _pEigVecs == NULL )
    error("[E:%s:%d %s] Invalid eigenvectors", __FILE__, __LINE__, __PRETTY_FUNCTION__ );

  pEigVecs = _pEigVecs;
  nsamples = (int32_t)pEigVecs->rows();
  ndims = (int32_t)pEigVecs->cols();

  pSVD = new Eigen::BDCSVD<Eigen::MatrixXd>(*_pEigVecs, Eigen::ComputeThinU | Eigen::ComputeThinV);

  // set dimensions;
  tol = _tol;
  maxLambda = _maxLambda;  
  
  hdr = NULL; wdr = NULL; iv = NULL;
  hwe0z = hwe1z = ibc0 = ibc1 = 0;

  pls = NULL; n_pls = 0; ploidies = NULL;
  ifs = new float[nsamples];
  betas = new float[ndims];
  pooled_af = -1;
  isaf_computed = false;

  gt = NULL;
  gq = NULL;
}

frequency_estimator::~frequency_estimator() {
  if ( ( pEigVecs ) && ( pSVD ) )
    delete pSVD;
  if (ifs != NULL )
    delete [] ifs;
  if ( gt != NULL )
    free(gt);
  if ( gq != NULL )
    free(gq);
}

frequency_estimator::frequency_estimator(Eigen::BDCSVD<Eigen::MatrixXd>* _pSVD, double _tol, double _maxLambda) :
  skipIf(false), skipInfo(false), siteOnly(false), gtError(0.005)   
{
  // set dimensions;
  pEigVecs = NULL;

  pSVD = _pSVD;
  nsamples = (int32_t)pSVD->matrixU().rows();
  ndims = (int32_t)pSVD->matrixU().cols();

  tol = _tol;
  maxLambda = _maxLambda;    
  
  hdr = NULL; wdr = NULL; iv = NULL;
  hwe0z = hwe1z = ibc0 = ibc1 = 0;

  pls = NULL; n_pls = 0; ploidies = NULL;
  ifs = new float[nsamples];
  betas = new float[ndims];
  pooled_af = -1;
  isaf_computed = false;  
}

bool frequency_estimator::set_hdr(bcf_hdr_t* _hdr, bcf_hdr_t* _wdr ) {
  if ( hdr != _hdr ) {
    hdr = _hdr;
    wdr = _wdr;
    char buffer[65535];
    if ( !skipInfo ) {
      if ( bcf_hdr_id2int(_hdr, BCF_DT_ID, "HWE_SLP_P" ) < 0 ) {
	sprintf(buffer,"##INFO=<ID=HWE_SLP_P,Number=1,Type=Float,Description=\"Z-score of HWE test with pooled allele frequencyes\">\n");
	bcf_hdr_append(wdr, buffer);
      }
      if ( bcf_hdr_id2int(_hdr, BCF_DT_ID, "IBC_P" ) < 0 ) {      
	sprintf(buffer,"##INFO=<ID=IBC_P,Number=1,Type=Float,Description=\"Inbreeding coefficient with pooled allele frequencies\">\n");
	bcf_hdr_append(wdr, buffer);
      }
      if ( bcf_hdr_id2int(_hdr, BCF_DT_ID, "HWE_SLP_I" ) < 0 ) {      
	sprintf(buffer,"##INFO=<ID=HWE_SLP_I,Number=1,Type=Float,Description=\"Z-score of HWE test with individual-sepcific allele frequencyes\">\n");
	bcf_hdr_append(wdr, buffer);
      }
      if ( bcf_hdr_id2int(_hdr, BCF_DT_ID, "IBC_I" ) < 0 ) {
	sprintf(buffer,"##INFO=<ID=IBC_I,Number=1,Type=Float,Description=\"Inbreeding coefficient with individual-sepcific allele frequencies\">\n");
	bcf_hdr_append(wdr, buffer);
      }
      if ( bcf_hdr_id2int(_hdr, BCF_DT_ID, "MAX_IF" ) < 0 ) {
	sprintf(buffer,"##INFO=<ID=MAX_IF,Number=1,Type=Float,Description=\"Maximum Individual-specific allele frequency\">\n");
	bcf_hdr_append(wdr, buffer);
      }
      if ( bcf_hdr_id2int(_hdr, BCF_DT_ID, "MIN_IF" ) < 0 ) {      
	sprintf(buffer,"##INFO=<ID=MIN_IF,Number=1,Type=Float,Description=\"Minimum Individual-specific allele frequency\">\n");
	bcf_hdr_append(wdr, buffer);
      }
      if ( bcf_hdr_id2int(_hdr, BCF_DT_ID, "BETA_IF" ) < 0 ) {      
	sprintf(buffer,"##INFO=<ID=BETA_IF,Number=%d,Type=Float,Description=\"Coefficients for intercept and each eigenvector to obtain ISAF\">\n", ndims);
	bcf_hdr_append(wdr, buffer);
      }
    }
    if ( ( !skipIf ) && ( !siteOnly ) ) {
      sprintf(buffer,"##FORMAT=<ID=IF,Number=1,Type=Float,Description=\"Individual-specific allele frequencies\">\n");
      bcf_hdr_append(wdr, buffer);
    }
    bcf_hdr_sync(wdr);
    return true;
  }
  else return false;
}

bool frequency_estimator::set_variant(bcf1_t* _iv, int8_t* _ploidies, int32_t* _pl) { //, std::vector<int32_t>* p_icols) {
  iv = _iv;
  ploidies = _ploidies;
  if ( iv->n_sample != nsamples )
    error("[E:%s:%d %s] nsamples %d != %d in the EigenVector", __FILE__, __LINE__, __PRETTY_FUNCTION__, iv->n_sample, nsamples);

  //error("%d %d %d",bcf_hdr_nsamples(hdr),iv->n_sample,nsamples);      
    
  // parse PL fields
  bcf_unpack(iv, BCF_UN_ALL);
    
  if ( _pl != NULL ) { pls = _pl; n_pls = 3*nsamples; }
  else {
    bool plfound = false;
    
    if ( field.empty() || ( field.compare("PL") == 0 ) ) {
      if ( bcf_get_format_int32(hdr, iv, "PL", &pls, &n_pls) < 0 ) {
	if ( field.compare("PL") == 0 ) 
	  error("[E:%s:%d %s] Cannot parse PL field", __FILE__, __LINE__, __PRETTY_FUNCTION__);
      }
      else
	plfound = true;
    }
    
    if ( ( (!plfound) && field.empty() ) || field.compare("GL") == 0 ){
      float* gls = NULL;
      int32_t n_gls = 0;
      if ( bcf_get_format_float(hdr, iv, "GL", &gls, &n_gls) < 0 ) {
	error("[E:%s:%d %s] Cannot parse GL field", __FILE__, __LINE__, __PRETTY_FUNCTION__);
      }
      else {
	if ( pls == NULL ) pls = new int32_t[n_gls];      
	for(int32_t i=0; i < nsamples; ++i) {
	  float maxgl = gls[3*i];
	  if ( gls[3*i+1] > maxgl ) maxgl = gls[3*i+1];
	  if ( gls[3*i+2] > maxgl ) maxgl = gls[3*i+2];
	  pls[3*i] = (int32_t)floor(-10*(gls[3*i]-maxgl)+0.5);
	  pls[3*i+1] = (int32_t)floor(-10*(gls[3*i+1]-maxgl)+0.5);
	  pls[3*i+2] = (int32_t)floor(-10*(gls[3*i+2]-maxgl)+0.5);
	  if ( pls[3*i] > 255 ) pls[3*i] = 255;
	  if ( pls[3*i+1] > 255 ) pls[3*i+1] = 255;
	  if ( pls[3*i+2] > 255 ) pls[3*i+2] = 255;
	}
	free(gls);
	n_pls = n_gls;
      }
    }
    else if ( field.compare("GT") == 0 ) {
      int32_t* gts = NULL;
      int32_t n_gts = 0;
      
      double tmp = gtError + gtError*gtError;
      int32_t errM[9] =
	{ 0, (int32_t)floor(-10*log10(gtError/(1-tmp))+0.5), (int32_t)floor(-10*log10(gtError*gtError/(1-tmp))+0.5),
	  (int32_t)floor(-10*log10(tmp/(2-tmp-tmp))+0.5), 0, 0,
	  0, 0, 0 };
      errM[5] = errM[3];
      errM[6] = errM[2];    
      errM[7] = errM[1];
      
      if ( bcf_get_genotypes(hdr, iv, &gts, &n_gts) < 0 ) {
	error("[E:%s:%d %s] Cannot parse GT field", __FILE__, __LINE__, __PRETTY_FUNCTION__);
      }
      else {
	int32_t max_ploidy = n_gts/nsamples;
	if ( max_ploidy != 2 )
	  error("[E:%s:%d %s] Multi-allelic (or Mono-allelic) variants found", __FILE__, __LINE__, __PRETTY_FUNCTION__);	
	if ( pls == NULL ) pls = new int32_t[nsamples*3];      
	for(int32_t i=0; i < nsamples; ++i) {
	  int32_t geno = bcf_gt_allele(gts[2*i])+bcf_gt_allele(gts[2*i+1]);	
	  if ( bcf_gt_is_missing(gts[2*i+1]) ) { // haploid or missing
	    if ( bcf_gt_is_missing(gts[2*i]) ) { // missing
	      pls[3*i] = pls[3*i+1] = pls[3*i+2] = 0;
	      continue;
	    }
	  else { // pretend to be homozygous for haploid
	    geno = bcf_gt_allele(gts[2*i]) + bcf_gt_allele(gts[2*i]);
	  }
	  }
	  pls[3*i] = errM[geno*3];
	  pls[3*i+1] = errM[geno*3+1];
	  pls[3*i+2] = errM[geno*3+2];
	}
	free(gts);
	n_pls = nsamples*3;
      }
    }
    else {
      if ( !plfound ) 
	error("[E:%s:%d %s] Cannot recognize the field [%s]", __FILE__, __LINE__, __PRETTY_FUNCTION__, field.c_str());
    }
  }
  /*
  else if ( bcf_get_format_int32(hdr, iv, "PL", &pls, &n_pls) < 0 ) {
    //float gls[nsamples+1];
    float* gls = NULL;
    int32_t n_gls = 0;
    if ( bcf_get_format_float(hdr, iv, "GL", &gls, &n_gls) < 0 ) {
      error("[E:%s:%d %s] Cannot parse PL or GL field", __FILE__, __LINE__, __PRETTY_FUNCTION__);
    }
    else {
      if ( pls == NULL ) pls = new int32_t[n_gls];      
      for(int32_t i=0; i < nsamples; ++i) {
	float maxgl = gls[3*i];
	if ( gls[3*i+1] > maxgl ) maxgl = gls[3*i+1];
	if ( gls[3*i+2] > maxgl ) maxgl = gls[3*i+2];
	pls[3*i] = (int32_t)floor(-10*(gls[3*i]-maxgl)+0.5);
	pls[3*i+1] = (int32_t)floor(-10*(gls[3*i+1]-maxgl)+0.5);
	pls[3*i+2] = (int32_t)floor(-10*(gls[3*i+2]-maxgl)+0.5);
	if ( pls[3*i] > 255 ) pls[3*i] = 255;
	if ( pls[3*i+1] > 255 ) pls[3*i+1] = 255;
	if ( pls[3*i+2] > 255 ) pls[3*i+2] = 255;
      }
      free(gls);      
    }
  }
  */
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
    U0 = pooled_af*(1-pooled_af)*(l0-2*l1+l2)/sum0;
    sumU0 += U0;
    sqU0 += (U0*U0);
    
    if ( use_isaf ) {
      estimate_isaf_em();
      
      sum1 = l0*(1.-ifs[i])*(1.-ifs[i]) + 2*l1*(1.-ifs[i])*ifs[i] + l2*ifs[i]*ifs[i] + 1e-100;
      ph1 = 2*l1*(1.-ifs[i])*ifs[i];
      obsH1 += (ph1/sum1);
      expH1 += (2*(1.-ifs[i])*ifs[i]);
      U1 = (1-ifs[i])*ifs[i]*(l0-2*l1+l2)/sum1;
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
	  //UD2(j,k) *= ( d2[k] / ( d2[k] + lambda ) );
	  UD2(j,k) *= ( d2[k] * d2[k] / ( d2[k] + lambda ) / ( d2[k] + lambda) );
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

    Eigen::VectorXd d2 = pSVD->singularValues();    
    Eigen::MatrixXd VD = pSVD->matrixV();
    for(int32_t j=0; j < ndims; ++j) {
      for(int32_t k=0; k < ndims; ++k) {
	VD(j,k) *= ( d2[k] / ( d2[k] + lambda ) / ( d2[k] + lambda) );	
      }
    }
    Eigen::VectorXd vBeta = VD * ( pSVD->matrixU().transpose() * y ) / 2.0;
    for(int32_t k=0; k < ndims; ++k)
      betas[k] = (float)vBeta[k];
    
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

bool frequency_estimator::update_gt_gq(bool update_gq) {
  if ( siteOnly ) return false;
  double gp = 0, gp_sum = 0, max_gp = 0;
  int32_t best_gt = 0;
  int32_t best_a1 = 0, best_a2 = 0;
  int32_t an = 0;
  int32_t acs[2] = {0,0};
  int32_t gcs[3] = {0,0,0};
  float afs[3];
  int32_t max_gq = 0;

  if ( gt == NULL ) 
    gt = (int32_t*) malloc(sizeof(int32_t)*2*nsamples);
  if ( ( update_gq ) && ( gq == NULL ) ) 
    gq = (int32_t*) malloc(sizeof(int32_t)*nsamples);    
  
  for(int32_t i=0; i < nsamples; ++i) {
    int32_t* pli = &pls[ i * 3 ];
    
    if ( ploidies[i] == 1 ) {
      max_gp = gp_sum = gp = ( phredConv.toProb((uint32_t)pli[0]) * (1.0 - ifs[i]) );
      best_gt = 0; best_a1 = 0; best_a2 = 0;
      gp = ( phredConv.toProb((uint32_t)pli[2]) * ifs[i] );
      gp_sum += gp;
      if ( max_gp < gp ) {
	max_gp = gp;
	best_gt = 2; best_a1 = 1; best_a2 = 1;
      }      
    }
    else if ( ploidies[i] == 2 ) {
      max_gp = gp_sum = gp = ( phredConv.toProb((uint32_t)pli[0]) * (1.0-ifs[i]) * (1.0-ifs[i]) );
      best_gt = 0; best_a1 = 0; best_a2 = 0;

      gp = phredConv.toProb((uint32_t)pli[1]) * 2.0 * ifs[i] * (1.0-ifs[i]);
      gp_sum += gp;
      if ( max_gp < gp ) { max_gp = gp; best_gt = 1; best_a1 = 0; best_a2 = 1; }

      gp = phredConv.toProb((uint32_t)pli[2]) * ifs[i] * ifs[i];
      gp_sum += gp;
      if ( max_gp < gp ) { max_gp = gp; best_gt = 2; best_a1 = 1; best_a2 = 1; }      
    }
    else if ( ploidies[i] == 0 ) {
      best_gt = 0;
      max_gp = 0;
      gp_sum = 1e-100;      
    }
    else
      error("[E:%s:%d %s] Unexpected ploidy %d", __FILE__, __LINE__, __PRETTY_FUNCTION__, (int32_t)ploidies[i]);

    if ( update_gq ) {
      double prob = 1.-max_gp/gp_sum;  // to calculate GQ
      if ( prob <= 3.162278e-26 )
	prob = 3.162278e-26;
      if ( prob > 1 )
	prob = 1;
    
      gq[i] = (int32_t)phredConv.err2Phred((double)prob);
      if ( ( best_gt > 0 ) && ( max_gq < gq[i] ) ) {
	max_gq = gq[i];
      }
    }
    
    gt[2*i]   = ((best_a1 + 1) << 1);
    gt[2*i+1] = ((best_a2 + 1) << 1);	    
    an += 2;             // still use diploid representation of chrX for now.
    ++acs[best_a1];
    ++acs[best_a2];
    ++gcs[best_gt];    
  }

  for(size_t i=0; i < 2; ++i) {
    afs[i] = acs[i]/(float)an;
  }

  //notice("Calling bcf_update_format_int32() with nsamples=%d",nsamples);

  bcf_update_format_int32(hdr, iv, "GT", gt, nsamples * 2);
  if ( update_gq )
    bcf_update_format_int32(hdr, iv, "GQ", gq, nsamples );	  

  iv->qual = (float)max_gq;

  bcf_update_info_int32(hdr, iv, "AC", &acs[1], 1);
  bcf_update_info_int32(hdr, iv, "AN", &an, 1);
  bcf_update_info_float(hdr, iv, "AF", &afs[1], 1);
  bcf_update_info_int32(hdr, iv, "GC", gcs, 3);
  bcf_update_info_int32(hdr, iv, "GN", &nsamples, 1);

  return true;
}
 
bool frequency_estimator::update_variant() {
  float hweslp0 = (float)((hwe0z > 0 ? -1 : 1) * log10( erfc(fabs(hwe0z)/sqrt(2.0)) + 1e-100 ));
  float hweslp1 = (float)((hwe1z > 0 ? -1 : 1) * log10( erfc(fabs(hwe1z)/sqrt(2.0)) + 1e-100 ));
  float max_if = 0, min_if = 1;
  for(int32_t j=0; j < nsamples; ++j) {
    if ( ifs[j] > max_if ) max_if = ifs[j];
    if ( ifs[j] < min_if ) min_if = ifs[j];
  }

  if ( siteOnly ) {
    //notice("foo");
    bcf_subset(hdr, iv, 0, 0);
    //notice("goo");    
  }  

  if ( !skipInfo ) {
    float hweaf = (float)pooled_af;
    bcf_update_info_float(wdr, iv, "HWEAF_P", &hweaf, 1);
    bcf_update_info_float(wdr, iv, "FIBC_P", &ibc0, 1);    
    bcf_update_info_float(wdr, iv, "HWE_SLP_P", &hweslp0, 1);
    bcf_update_info_float(wdr, iv, "FIBC_I", &ibc1, 1);    
    bcf_update_info_float(wdr, iv, "HWE_SLP_I", &hweslp1, 1);
    bcf_update_info_float(wdr, iv, "MAX_IF", &max_if, 1);
    bcf_update_info_float(wdr, iv, "MIN_IF", &min_if, 1);
    bcf_update_info_float(wdr, iv, "BETA_IF", betas, ndims);
  }
  if ( ( !skipIf ) && ( !siteOnly ) ) {
    bcf_update_format_float(wdr, iv, "IF", ifs, nsamples);
  }
  
  return true;
}
