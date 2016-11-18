#include "cramore.h"

int32_t cmdScMultinomEM(int32_t argc, char** argv) {
  std::string inMatrix;
  std::string outPrefix;
  double alpha = 1.0;       // pseudo-count per cell
  double thresDiff = 1e-10; // threshold to stop EM iteration
  int32_t maxIter = 100;    // maximum number of EM iteration
  int32_t nClust = 0;       // Number of clusters required
  int32_t nRestarts = 1;    // Number of restarts to pick the best model
  int32_t seed = 0;         // random seed
  int32_t nCollapseGenes = 0; // collapse genes into a specific number
  double fracSubsample = 1;   // fraction of samples to thin the data

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Required Options", NULL)
    LONG_STRING_PARAM("in",&inMatrix, "Input matrix int the format of R-compatible text matrix (can be gzipped)")
    LONG_STRING_PARAM("out",&outPrefix, "Output file prefix")
    LONG_INT_PARAM("k",&nClust, "Number of clusters")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_DOUBLE_PARAM("alpha",&alpha, "Pseudo-count per cell")    
    LONG_DOUBLE_PARAM("thres",&thresDiff, "Threshold of LLK difference to terminate the EM iteration")
    LONG_INT_PARAM("restarts",&nRestarts, "Number of restarts to pick the best model")
    LONG_INT_PARAM("max-iter",&maxIter, "Number of maximum E-M iterations")
    LONG_INT_PARAM("collapse-genes",&nCollapseGenes,"Number of genes to be collapsed into to reduce parameter space")
    LONG_INT_PARAM("seed",&seed, "Seed for random number generator (default uses clock)")
    LONG_DOUBLE_PARAM("frac-subsample",&fracSubsample, "Fraction of samples to thin the data")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  if ( nClust == 0 ) {
    error("[E:%s:%d %s] --k is a required parameter",__FILE__,__LINE__,__FUNCTION__);
  }

  if ( outPrefix.empty() ) {
    error("[E:%s:%d %s] --out is a required parameter",__FILE__,__LINE__,__FUNCTION__);
  }

  if ( inMatrix.empty() ) {
    error("[E:%s:%d %s] --in is a required parameter",__FILE__,__LINE__,__FUNCTION__);
  }

  htsFile* wf = hts_open((outPrefix+".pis").c_str(),"w");
  if ( wf == NULL )
    error("[E:%s:%d %s] Cannot open file %s for writing",__FILE__,__LINE__,__FUNCTION__, (outPrefix+".pis").c_str());

  htsFile* hp = hts_open(inMatrix.c_str(), "r");
  if ( hp == NULL )
    error("[E:%s:%d %s] Cannot open file %s for reading",__FILE__,__LINE__,__FUNCTION__,inMatrix.c_str());

  kstring_t str = {0,0,0};  
  int32_t lstr = 0;

  // read and parse header columns
  lstr = hts_getline(hp, KS_SEP_LINE, &str);
  if ( lstr < 0 )
    error("[E:%s:%d %s] Cannot find header line from %s",__FILE__,__LINE__,__FUNCTION__,inMatrix.c_str());

  
  int32_t nfields = 0;
  int32_t* fields = NULL;  
  fields = ksplit(&str, 0, &nfields);

  std::vector<std::string> hdrs;
  for(int32_t i=0; i < nfields; ++i) {
    hdrs.push_back(std::string(&str.s[fields[i]]));
  }

  std::vector<int32_t*> R;
  std::vector<std::string> genes;
  std::vector<int64_t> rowSums;
  int64_t* colSums = NULL;

  notice("%d columns found in the header",(int32_t)hdrs.size());

  // read and parse the matrix elements
  int64_t nZero = 0;
  int64_t nSum = 0;
  int32_t nEmptyRows = 0;
  while( ( lstr = hts_getline(hp, KS_SEP_LINE, &str) ) >= 0 ) {
    fields = ksplit(&str, 0, &nfields);

    if ( R.empty() ) {
      notice("%d columns found in the second line",nfields);      
      if ( nfields == (int32_t)hdrs.size() ) { // remove one from the header
	hdrs.erase(hdrs.begin());
	notice("Ignoring the first column in the header");	
      }
      colSums = (int64_t*)calloc((int32_t)hdrs.size(),sizeof(int64_t));
    }
    else {
      if ( ( nfields != (int32_t)hdrs.size() + 1 ) && ( nfields != (int32_t)hdrs.size() + 0 ) )
	error("[E:%s:%d %s] Inconsistent number of headers. Expected %d but observed %d",__FILE__,__LINE__,__FUNCTION__,(int32_t)hdrs.size()+1, nfields);
    }

    int32_t* cnts = (int32_t*)malloc(sizeof(int32_t)*(nfields-1));
    int64_t rowSum = 0;
    for(int32_t i=1; i < nfields; ++i) {
      cnts[i-1] = atoi(&str.s[fields[i]]);

      if ( fracSubsample < 1 ) {
	int32_t tot = cnts[i-1];
	int32_t sampled = 0;
	for(int32_t j=0; j < tot; ++j) {
	  if ( (rand()+0.5) / (RAND_MAX+1.) < fracSubsample )
	    ++sampled;
	}
	cnts[i-1] = sampled;
      }
      
      if ( cnts[i-1] == 0 ) ++nZero;
      else {
	rowSum += cnts[i-1];
	colSums[i-1] += cnts[i-1];
      }
    }
    
    if ( rowSum == 0 ) {
      free(cnts);
      ++nEmptyRows;
    }
    else {
      genes.push_back(std::string(&str.s[fields[0]]));      
      R.push_back(cnts);
      rowSums.push_back(rowSum);
      nSum += rowSum;
    }
  }
  hts_close(hp);


  int32_t nRow = (int32_t)genes.size();
  int32_t nCol = (int32_t)hdrs.size();
  int64_t nCell = (int64_t)nRow * (int64_t)nCol;

  notice("Loaded a matrix with %d rows and %d columns after ignoring %d empty rows. Sparsity is %.5lg. Average of non-empty cells is %.5lg", nRow, nCol, nEmptyRows, (double)nZero/(double)(nCell+nEmptyRows*nCol), (double)nSum/(double)(nCell+nEmptyRows*nCol-nZero));

  if ( nCollapseGenes > 0 ) {
    std::vector< std::vector<int32_t> > group2Gene( nCollapseGenes );
    std::vector< int32_t > gene2Group( nRow, 0 );

    for(int32_t i=0; i < nRow; ++i) {
      int32_t g = (rand() % nCollapseGenes);
      group2Gene[g].push_back(i);
    }

    nEmptyRows = 0;
		
    for(int32_t i=nCollapseGenes-1; i >= 0; --i) {
      if ( group2Gene[i].empty() ) {
	++nEmptyRows;
	group2Gene.erase(group2Gene.begin() + i);
      }
      else {
	for(int32_t j=0; j < (int32_t)group2Gene[i].size(); ++j) {
	  gene2Group[group2Gene[i][j]] = i;
	}
      }
    }

    std::vector<std::string> newGenes(nCollapseGenes-nEmptyRows);
    std::vector<int64_t> newRowSums(nCollapseGenes-nEmptyRows, 0);
    std::vector<int32_t*> newR(nCollapseGenes-nEmptyRows, NULL);

    for(int32_t i=0; i < nRow; ++i) {
      int32_t g = gene2Group[i];
      if ( newGenes[g].empty() ) {
	newGenes[g] = genes[i];
	newR[g] = (int32_t*) calloc(sizeof(int32_t), nCol);
      }
      else {
	newGenes[g] += ",";
	newGenes[g] += genes[i];
      }
      newRowSums[g] += rowSums[i];
      for(int32_t j=0; j < nCol; ++j) {
	newR[g][j] += R[i][j];
      }
      free(R[i]);
    }

    genes = newGenes;
    rowSums = newRowSums;
    R = newR;
    nRow = nCollapseGenes-nEmptyRows;
    nCell = (int64_t)nRow * (int64_t)nCol;    

    notice("Collapsed the matrix with %d rows and %d columns after ignoring %d additional empty rows created during the random collpaing procedure", nRow, nCol, nEmptyRows);
  }

  // calculate the global proportion matrix
  double* p0 = new double[nRow];
  for(int32_t i=0; i < nRow; ++i) {
    p0[i] = (double)rowSums[i]/(double)nSum;
  }

  // create multiple copies of parameters for simultaneous EM
  int32_t nCxR = nClust * nRestarts;
  double* pis = (double*)calloc(nCxR,sizeof(double));
  double* Ps = (double*)calloc(nCxR * nRow,sizeof(double));
  double* Zs = (double*)calloc(nCxR * nCol,sizeof(double));
  double* llks = new double[nRestarts];
  double* llk0s = new double[nRestarts];

  if ( seed == 0 )
    srand(time(NULL));
  else
    srand(seed);

  // randomize class assignments
  for(int32_t c=0; c < nCol; ++c) {
    double* z = &Zs[c * nCxR];
    for(int32_t r=0; r < nRestarts; ++r) {
      z[(int32_t)(floor((rand()+0.5)/(RAND_MAX+1.)*nClust)) * nRestarts + r] = 1.;
    }
  }

  //notice("foo");

  // run EM iteration
  for(int32_t iter=0; iter < maxIter; ++iter) {
    // M-step for pi
    for(int32_t c=0; c < nCol; ++c) {
      double* z = &Zs[c * nCxR];    
      for(int32_t k=0; k < nCxR; ++k) {
	pis[k] += z[k];  // pi_k = \sum_c Pr(z_c = k)
      }
    }

    //for(int32_t k=0; k < nCxR; ++k) {
    //  notice("k=%d\tu_pi=%lg",k,pis[k]);
    //}    

    // normalize pi_k, and take a log
    for(int32_t r=0; r < nRestarts; ++r) {
      double sum = 0;
      for(int32_t k=0; k < nClust; ++k) {
	sum += pis[k*nRestarts+r];
      }
      for(int32_t k=0; k < nClust; ++k) {
	pis[k*nRestarts+r] = log( pis[k*nRestarts+r]/sum );
      }	
    }

    //for(int32_t k=0; k < nCxR; ++k) {
    //  notice("k=%d\tpi=%lg",k,exp(pis[k]));
    //}

    //notice("bar");    
    
    // M-step for P (without normalization)
    for(int32_t g=0; g < nRow; ++g) {
      double* p = &Ps[g * nCxR];
      
      for(int32_t k=0; k < nCxR; ++k)
	p[k] = 0;
      
      for(int32_t c=0; c < nCol; ++c) {
	double* z = &Zs[c * nCxR];
	double r = R[g][c] + p0[g]*alpha;
	//double r = R[g][c] + alpha/nRow;
	for(int32_t k=0; k < nCxR; ++k) {
	  //double t = z[k] * r;
	  p[k] += (z[k] * r); // not normalized   \Pr(x_g|z_c=k) \propt \sum_c R_gc Pr(z_c=k)
	}
      }
    }

    //notice("goo");

    // normalize P
    for(int32_t k=0; k < nCxR; ++k) {
      double sumP = 0;
      for(int32_t g=0; g < nRow; ++g) {
	sumP += Ps[g*nCxR + k];
      }
      
      for(int32_t g=0; g < nRow; ++g) {
	Ps[g*nCxR + k] /= sumP;
      }
    }
    

    // transform p into logp
    for(int32_t g=0; g < nRow; ++g) {
      double* p = &(Ps[g * nCxR]);
      for(int32_t k=0; k < nCxR; ++k) {
	p[k] = log(p[k]); // + pis[k];  // pi*P in log-scale
      }
    }


    for(int32_t c=0; c < nCol; ++c) {
      double* z = &(Zs[c*nCxR]);      
      for(int32_t k=0; k < nCxR; ++k) {
	z[k] = pis[k];
      }      
    }
    
    
    // E-step for Z : t(R) %*% logP
    for(int32_t g=0; g < nRow; ++g) {
      double* p = &(Ps[g * nCxR]);      
      for(int32_t c=0; c < nCol; ++c) {
	int32_t r = R[g][c];
	double* z = &(Zs[c*nCxR]);
	for(int32_t k=0; k < nCxR; ++k) {
	  z[k] += (r*p[k]);  // \log Pr(z_c = k | x) \propt \sum_g [ \log \Pr(R_gc|z_c=k) ] + \pi_k
	}
      }
    }

    for(int32_t r=0; r < nRestarts; ++r) {
      if ( iter == 0 )
	llk0s[r] = -1e300;
      else
	llk0s[r] = llks[r];
      
      llks[r] = 0;      
    }

    for(int32_t c=0; c < nCol; ++c) {
      double* z = &(Zs[c*nCxR]);      
      for(int32_t r=0; r < nRestarts; ++r) {
	double maxZ = z[r];
	for(int32_t k=1; k < nClust; ++k) {
	  if ( maxZ < z[k*nRestarts + r])
	    maxZ = z[k*nRestarts + r];
	}

	double sumZ = 0;
	for(int32_t k=0; k < nClust; ++k) {
	  double zdiff = z[k*nRestarts+r] - maxZ;	  
	  //sumZ += (z[k*nRestarts+r] = ((zdiff < -100) ? 3.7e-44 : exp(zdiff)));
	  sumZ += (z[k*nRestarts+r] = exp(zdiff));
	}

	//for(int32_t k=0; k < nClust; ++k) {
	//  z[k*nRestarts+r] /= sumZ;
	//}
	
	llks[r] += (maxZ + log(sumZ));
	//notice("%d\t%lg\t%lg\t%lg\t%lg\t%lg",c,llks[r],maxZ,sumZ,z[0],z[1]);
	if ( isnan(llks[r]) )
	  abort();
      }
    }
    double maxDiff = -1e300;
    for(int32_t r=0; r < nRestarts; ++r) {
      notice("Iter:%d\tThread %d:\tLLK=%.5lf\tDiff=%.5lg", iter, r, llks[r], llks[r]-llk0s[r]); //, exp(pis[0]), exp(pis[nRestarts]));
      if ( maxDiff < llks[r]-llk0s[r] ) {
	maxDiff = llks[r]-llk0s[r];
      }
    }

    if ( maxDiff < thresDiff ) {
      notice("All LLK differences are less than %.5lg < %.5lg", maxDiff, thresDiff);
      break;
    }

    if ( iter + 1 == maxIter )
      notice("Reached maximum iteration %d",maxIter);      
  }

  int32_t iMin = 0;
  for(int32_t r=1; r < nRestarts; ++r) {
    if ( llks[iMin] < llks[r] )
      iMin = r;
  }

  // transform P to linear scale
  for(int32_t k=0; k < nClust; ++k) {
    int32_t kr = k*nRestarts+iMin;
    double maxP = Ps[kr];
    for(int32_t g=0; g < nRow; ++g) {
      if ( maxP < Ps[g*nCxR + kr] )
	maxP = Ps[g*nCxR + kr];
    }

    double sumP = 0;
    for(int32_t g=0; g < nRow; ++g) {
      sumP += (Ps[g*nCxR + kr] = exp(Ps[g*nCxR + kr] - maxP));
    }

    for(int32_t g=0; g < nRow; ++g) {
      Ps[g*nCxR + kr] /= sumP;
    }
  }

  for(int32_t k=0; k < nClust; ++k) 
    hprintf(wf, "%g\n",exp(pis[k*nRestarts+iMin]));
  hts_close(wf);

  wf = hts_open((outPrefix+".Ps").c_str(),"w");
  for(int32_t g=0; g < nRow; ++g) {
    hprintf(wf, "%s",genes[g].c_str());    
    for(int32_t k=0; k < nClust; ++k) {
      hprintf(wf, "\t%.5lg",Ps[g*nCxR + k*nRestarts + iMin]); 
    }
    hprintf(wf, "\n");
  }
  hts_close(wf);

  wf = hts_open((outPrefix+".Zs").c_str(),"w");
  for(int32_t c=0; c < nCol; ++c) {
    double sumZ = 0;
    for(int32_t k=0; k < nClust; ++k)
      sumZ += Zs[c*nCxR + k*nRestarts + iMin];

    int32_t iBest = 0;
    for(int32_t k=1; k < nClust; ++k) {
      if ( Zs[c*nCxR + iBest*nRestarts + iMin] < Zs[c*nCxR + k*nRestarts + iMin] )
	iBest = k;
    }

    hprintf(wf, "%s\t%d\t%d",hdrs[c].c_str(), colSums[c], iBest+1);
    for(int32_t k=0; k < nClust; ++k)
      hprintf(wf, "\t%.5lg",Zs[c*nCxR + k*nRestarts + iMin]/sumZ);
    hprintf(wf, "\n");
  }
  hts_close(wf);    
  
  // free up the memories
  for(int32_t i=0; i < nRow; ++i) {
    free(R[i]);
  }
  delete[] llks;
  delete[] p0;
  free(pis);
  free(Zs);
  free(Ps);
  free(colSums);
  
  return 0;
}

