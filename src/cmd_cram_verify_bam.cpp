#include "cramore.h"
#include "bcf_filtered_reader.h"
#include "sam_filtered_reader.h"
#include "var_dict.h"
#include "contam_estimator.h"

int32_t cmdCramVerifyBam(int32_t argc, char** argv) {
  SAMFilteredReader sr;
  int32_t capBQ = 40;
  int32_t minBQ = 13;
  int32_t minDP = 2;
  int32_t capDP = 50;
  int32_t minTD = 2;
  sr.filt.exclude_flag = 0x0f04;
  sr.filt.minMQ = 20;
  std::string tagRG = "RG";
  
  BCFFilteredReader vr;
  std::string field("GT");
  double genoError = 0.01;
  vr.vfilt.minMAC = 1;
  vr.vfilt.minCallRate = 0.5;
  vr.vfilt.maxAlleles = 2;

  std::string svdPrefix;
  int32_t numPC = 2;
  
  int32_t seed = 0;  

  std::string outPrefix;
  bool ignoreRG = false;
  std::string fixPC;
  double fixAlpha = -1;
  bool withinAncestry = false;

  bool refBias = false;
  bool writePileup = false;
  bool skipDup = false;

  std::vector<std::string> smIDs;

  vr.verbose = 1000000;
  sr.verbose = 10000000;  

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Options for input SAM/BAM/CRAM", NULL)
    LONG_STRING_PARAM("sam",&sr.sam_file_name, "Input SAM/BAM/CRAM file. Must be sorted by coordinates and indexed")
    LONG_INT_PARAM("cap-BQ", &capBQ, "Maximum base quality (higher BQ will be capped)")
    LONG_INT_PARAM("min-BQ", &minBQ, "Minimum base quality to consider (lower BQ will be skipped)")
    LONG_INT_PARAM("min-MQ", &sr.filt.minMQ, "Minimum mapping quality to consider (lower MQ will be ignored)")
    LONG_INT_PARAM("min-TD", &minTD, "Minimum distance to the tail (lower will be ignored)")
    LONG_INT_PARAM("excl-flag", &sr.filt.exclude_flag, "SAM/BAM FLAGs to be excluded")
    LONG_STRING_PARAM("tag-rg", &tagRG, "Tag representing readgroup")
    LONG_PARAM("ignore-rg",&ignoreRG, "Do not group by readgroups")
    LONG_PARAM("skip-dup",&skipDup, "Skip potential duplicate aggressively based on start position")
    LONG_INT_PARAM("min-DP", &minDP, "Minimum depth of reads with coverage to include the site")
    LONG_INT_PARAM("cap-DP", &capDP, "Maximum depth of reads with coverage to cap the depth")    
    

    LONG_PARAM_GROUP("Options for input VCF/BCF", NULL)
    LONG_STRING_PARAM("vcf",&vr.bcf_file_name, "Input VCF/BCF file, containing the individual genotypes (GT), posterior probability (GP), or genotype likelihood (PL)")
    LONG_STRING_PARAM("field",&field,"FORMAT field to extract the genotype, likelihood, or posterior from")
    LONG_DOUBLE_PARAM("geno-error",&genoError,"Genotype error rate (must be used with --field GT)")
    LONG_INT_PARAM("min-mac",&vr.vfilt.minMAC, "Minimum minor allele frequency")
    LONG_DOUBLE_PARAM("min-callrate",&vr.vfilt.minCallRate, "Minimum call rate")    
    LONG_MULTI_STRING_PARAM("sm",&smIDs, "List of sample IDs to compare to (default: use all)")
    LONG_STRING_PARAM("sm-list",&vr.sample_id_list, "File containing the list of sample IDs to compare")

    LONG_PARAM_GROUP("Options for input SVD files", NULL)
    LONG_STRING_PARAM("svd", &svdPrefix, "Prefix of SVD (.fUD, .V) files")
    LONG_INT_PARAM("num-PC", &numPC, "Number of PCs to use")    

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out",&outPrefix,"Output file prefix")
    LONG_STRING_PARAM("fix-pc", &fixPC, "Fix principal component coordinates in [x,y] format")
    LONG_DOUBLE_PARAM("fix-alpha", &fixAlpha, "Fix the contamination value")
    LONG_PARAM("within-ancestry", &withinAncestry, "Make the ancestries of intended and contaminating samples the same")
    LONG_PARAM("ref-bias", &refBias, "Estimate the reference bias parameters together")
    LONG_PARAM("write-pileup", &writePileup, "Writeup pileup information for each allele")

    LONG_PARAM_GROUP("Other toptions", NULL)    
    LONG_INT_PARAM("sam-verbose",&sr.verbose, "Verbose message frequency for SAM/BAM/CRAM")
    LONG_INT_PARAM("vcf-verbose",&vr.verbose, "Verbose message frequency for VCF/BCF")
    LONG_INT_PARAM("seed",&seed, "Random seed")    
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  if ( outPrefix.empty() || sr.sam_file_name.empty() || svdPrefix.empty() )
    error("Missing required options --out, --svd or --sam");

  if ( (int32_t)svdPrefix.empty() + (int32_t)vr.bcf_file_name.empty() != 1 )
    error("Only one of --svd or --vcf options are required");

  std::vector<std::string> refIDs;
  std::vector< std::vector<double> > matv;

  /*
  if ( !vr.bcf_file_name.empty() ) {
    for(int32_t i=0; i < (int32_t)smIDs.size(); ++i) {
      vr.add_specified_sample(smIDs[i].c_str());
    }

    vr.unlimited_buffer = true;
    vr.vfilt.maxAlleles = 2;
    vr.init_params();
  }
  */
  
  tsv_reader tsv_svd_u((svdPrefix+".fUD.gz").c_str());
  int32_t ncols = tsv_svd_u.read_line();
  if ( ncols < numPC + 2 )
    error("[E:%s:%d %s] observed %d < %d+2 columns",__FILE__,__LINE__,__PRETTY_FUNCTION__, ncols, numPC);

  var_dict<vb_plp> var2u;
  while( ( ncols = tsv_svd_u.read_line() ) > 0 ) {
    if ( ncols < numPC + 2 )
	error("[E:%s:%d %s] observed %d < %d+2 columns in the file",__FILE__,__LINE__,__PRETTY_FUNCTION__, ncols, numPC);
      
    double* v = new double[numPC+1];
    //std::vector<double>& v = var2u[str_field_at(0)];
    //v.resize(numPC+1);
    for(int32_t i=0; i <= numPC; ++i) {
      v[i] = tsv_svd_u.double_field_at(i+1);
    }
    var2u[tsv_svd_u.str_field_at(0)].ud = v;
  }
  
  tsv_reader tsv_svd_v((svdPrefix+".V.gz").c_str());
  ncols = tsv_svd_v.read_line();
  if ( ncols < numPC + 1 )
    error("[E:%s:%d %s] observed %d < %d+1 columns",__FILE__,__LINE__,__PRETTY_FUNCTION__, ncols, numPC);
  
  while( ( ncols = tsv_svd_v.read_line() ) > 0 ) {
    if ( ncols < numPC + 1 )
      error("[E:%s:%d %s] observed %d < %d+1 columns in the file",__FILE__,__LINE__,__PRETTY_FUNCTION__, ncols, numPC);
    
    matv.resize( matv.size() + 1 );
    std::vector<double>& v = matv.back();
    v.resize(numPC);
    refIDs.push_back(tsv_svd_v.str_field_at(0));
    for(int32_t i=0; i < numPC; ++i) {
      v[i] = tsv_svd_v.double_field_at(i+1);
    }
  }

  // calculate the mean and variance of PCs
  std::vector<double> sumv(numPC,0);
  std::vector<double> ssqv(numPC,0);
  std::vector<double> muv(numPC,0);
  std::vector<double> sdv(numPC,0);    
  for(int32_t i=0; i < (int32_t)matv.size(); ++i) {
    for(int32_t j=0; j < numPC; ++j) {
      sumv[j] += matv[i][j];
      ssqv[j] += (matv[i][j] * matv[i][j]);
    }
  }

  for(int32_t i=0; i < numPC; ++i) {
    muv[i] = sumv[i] / (int32_t)matv.size();
    sdv[i] = sqrt(ssqv[i]/(int32_t)matv.size() - muv[i]*muv[i]);
  }

  double minMAF = 0.25 / (double) matv.size();
    
  sr.set_buffer_size(1);
  sr.init_params();

  int32_t n_warning_no_gtag = 0;
  
  if ( outPrefix.empty() )
    error("[E:%s:%d %s] --out parameter is missing",__FILE__,__LINE__,__PRETTY_FUNCTION__);

  htsFile* wf = hts_open(outPrefix.c_str(), "w");
  if ( wf == NULL )
    error("Cannot open %s for writing", outPrefix.c_str());  

  char gtag[2] = {0,0};

  if ( ignoreRG || tagRG.empty() ) { // do nothing
  }
  else if ( tagRG.size() == 2 ) {
    gtag[0] = tagRG.at(0);
    gtag[1] = tagRG.at(1);    
  }
  else {
    error("[E:%s:%d %s] Cannot recognize group tag %s. It is suppose to be a length 2 string",__FILE__,__LINE__,__FUNCTION__,tagRG.c_str());
  }

  char base, qual;
  uint8_t allele, bq;
  int32_t rpos;
  kstring_t readseq = {0,0,0};
  kstring_t readqual = {0,0,0};

  //int32_t nReadsMultiSNPs = 0, nReadsSkipBCD = 0, nReadsPass = 0, nReadsRedundant = 0, nReadsN = 0, nReadsLQ = 0, nReadsTMP = 0;

  std::string it_chr;
  var_dict<vb_plp>::var_it_t it, it_l, it_r, it_beg, it_end;
  int32_t nreads = 0, nskips = 0, ncaps = 0;
  int32_t prevpos = 0;

  while( sr.read() ) { // read SAM file
    bam1_t* b = sr.cursor();    
    int32_t endpos = bam_endpos(b);
    int32_t begpos = b->core.pos;
    const char* chrom = bam_get_chrom(sr.hdr, b);

    if ( skipDup ) {
      if ( prevpos == begpos ) continue;
      prevpos = begpos;
    }
    
    if ( it_chr != chrom ) {
      it_chr = chrom;
      it_l = var2u.lower_bound(chrom, begpos);
      it_r = var2u.lower_bound(chrom, endpos);
      it_end = var2u.get_chr_end(chrom);
      it_beg = var2u.get_chr_beg(chrom);      
    }
    else {
      while( ( it_l != it_end ) && ( it_l->first.pos1 < begpos ) ) ++it_l;
      while( ( it_r != it_end ) && ( it_r->first.pos1 < endpos ) ) ++it_r;

      if ( it_l != it_r ) {
	//const char* rg = ".";
	if ( !ignoreRG && !tagRG.empty() ) {
	  uint8_t *tag = (*gtag) ? (uint8_t*) bam_aux_get(b, gtag) : NULL;
	  if ( ( tag != NULL ) && ( *tag == 'Z' ) ) {
	    //rg = bam_aux2Z(tag);
	  }
	  else {
	    if ( n_warning_no_gtag < 10 ) {
	      notice("WARNING: Cannot find RG tag %s from %d-th read %s at %s:%d-%d. Treating all of them as a single group", tagRG.c_str(), sr.n_read, bam_get_qname(b), chrom, begpos, endpos);
	    }
	    else if ( n_warning_no_gtag == 10 ) {
	      notice("WARNING: Suppressing 10+ missing RG tag warnings...");
	    }
	    ++n_warning_no_gtag;	
	  }
	}

	for(it = it_l; it != it_r; ++it) {
	  bam_get_base_and_qual_and_read_and_qual(b, (int32_t)it->first.pos1, base, qual, rpos, &readseq, &readqual);
	  if ( rpos == BAM_READ_INDEX_NA ) {
	    //if ( rand() % 1000 == 0 ) 
	    //notice("Cannot find any informative read between %s:%d-%d at %s:%d", bam_get_chrom(sr.hdr, b), b->core.pos+1, bam_endpos(b), bcf_hdr_id2name(vr.cdr.hdr, scl.snps[i].rid), scl.snps[i].pos);
	    continue;
	  }
	  if ( base == 'N' ) continue;
	  if ( it->first.alts.size() > 1 ) continue;

	  //++nv_valid;

	  if ( qual-33 < minBQ ) { continue; }
	  if ( rpos < minTD-1 ) { continue; }
	  if ( rpos + minTD > b->core.l_qseq ) { continue; }

	  allele = ( ( it->first.ref[0] == base ) ? 0 : ( ( it->first.alts[0][0] == base ) ? 1 : 2 ) );
	  bq = qual-33 > capBQ ? capBQ : qual-33;

	  if ( (int32_t)it->second.depth() < capDP ) {
	    it->second.add_al_bq(allele, bq);
	    //it->second.als.push_back(allele);
	    //it->second.bqs.push_back(bq);
	  }
	  else ++ncaps;
	}
	++nreads;

	if ( nreads % vr.verbose == 0 )
	  notice("Processed %d reads that overlap with sites in SVD", nreads);
      }
    }
    //if ( nv_valid > 0 ) ++scl.cell_totl_reads[ibcd];    
  }

  sr.close();

  if ( n_warning_no_gtag > 10 ) 
    notice("WARNING: Suppressed a total of %d RG warnings...", n_warning_no_gtag);

  // create a depth distribution
  std::vector<int32_t> dpCnts(capDP+1,0);
  var2u.get_depth_distribution(dpCnts);
  
  //notice("Finished reading %d reads across %d markers", nreads, var2u.nvar);

  contam_estimator est(numPC, minMAF);
  est.muPCs = muv;
  est.sdPCs = sdv;

  var2u.vectorize(est.plps);      // assign input data

  notice("Removing markers with less than %d reads..", minDP);

  int32_t v_rmv, v_kep, r_rmv, r_kep; 
  est.remove_low_depth(minDP, &v_rmv, &v_kep, &r_rmv, &r_kep); // clean the data

  notice("Summary of reads and variant sites");
  notice("Removed  reads    due to low depth : %d", r_rmv);
  notice("Removed  variants due to low depth : %d", v_rmv);
  notice("Retained reads    with sufficient depth : %d", r_kep);
  notice("Retained variants with sufficient depth : %d", v_kep);
  notice("Skipped  reads    above capped depth    : %d", ncaps);
  notice("Other skipped  reads                    : %d", nskips);
  
  //notice("Removed %d reads across %d markers and kept %d reads across %d markers above depth %d", nrm, (int32_t)est.plps.size(), minDP);  

  notice("Obtaining maximum likelihood estimates..", minDP);  

  double alpha;
  double* pc1 = new double[numPC];
  double* pc2 = new double[numPC];

  est.fixAlpha = fixAlpha;
  if ( fixPC.empty() ) est.fixPC = NULL;
  else {
    int32_t i;
    size_t j = 0, k = 0;

    est.fixPC = new double[numPC];
    for(i=0; i < numPC; ++i) {
      if ( ( k = fixPC.find(',',j) ) == std::string::npos ) {
	if ( i != numPC -1 )
	  error("Cannot parse %s as %d dimensional vector. Expecting more separator ','", fixPC.c_str(), numPC);
	k = (int32_t)fixPC.size();
      }
      else if ( i == numPC - 1 )
	error("Cannot parse %s as %d dimensional vector. More separators ',' observed as expected", fixPC.c_str(), numPC);
      //notice("%d %d %d %s %lf %s", i, j, k, fixPC.c_str(), atof("0,0"), fixPC.substr(j,k-j).c_str());
      est.fixPC[i] = atof(fixPC.substr(j,k-j).c_str());
      j = k+1;
      //notice("%d %d %d %lf", i, j, k, est.fixPC[i]);
    }
  }

  // OUTPUT FORMAT
  //
  // MODEL      KEY    VALUE
  // SAME_ANCESTRY   ALPHA  0.05
  // DIFF_ANCESTRY
  // FINAL_MODEL
  // BETWEEN  ALPHA  
  hprintf(wf, "ANCESTRY_MODEL\tKEY\tVALUE\n");
  // perform withinAncestry estimation first
  est.optimizeLLK(&alpha, pc1, pc2, true);

  double llk0 = est.llk0;
  double llk1 = est.llk1;

  int64_t sumCnts = 0, sumDepths = 0;
  for(int32_t i=0; i < (int32_t)dpCnts.size(); ++i) {
    sumCnts += dpCnts[i];
    sumDepths += (dpCnts[i] * i);
  }

  for(int32_t i=0; i < (int32_t)dpCnts.size(); ++i)  
    hprintf(wf, "DEPTH_COUNT\t%d\t%d\n", i, dpCnts[i]);
  
  for(int32_t i=0; i < (int32_t)dpCnts.size(); ++i)  
    hprintf(wf, "DEPTH_FRAC\t%d\t%.6lf\n", i, (double)dpCnts[i]/(double)(sumCnts+1e-10));

  hprintf(wf, "NOT_APPLICABLE\tMEAN_DEPTH\t%.6lf\n", (double)sumDepths/(double)(sumCnts+1e-10));

  hprintf(wf, "NOT_APPLICABLE\tREADS_KEPT\t%d\n",    r_kep);
  hprintf(wf, "NOT_APPLICABLE\tVARIANTS_KEPT\t%d\n", v_kep);
  hprintf(wf, "NOT_APPLICABLE\tREADS_LOW_DP\t%d\n",  r_rmv);
  hprintf(wf, "NOT_APPLICABLE\tVARIANTS_LOW_DP\t%d\n",v_rmv);
  hprintf(wf, "NOT_APPLICABLE\tREADS_CAPPED\t%d\n",  ncaps);
  hprintf(wf, "NOT_APPLICABLE\tREADS_SKIPPED\t%d\n", nskips);  
  
  hprintf(wf, "SAME_ANCESTRY\tLOG_LIKELIHOOD_H0\t%.5lg\n",est.llk0);
  hprintf(wf, "SAME_ANCESTRY\tLOG_LIKELIHOOD_H1\t%.5lg\n",est.llk1);
  hprintf(wf, "SAME_ANCESTRY\tCHISQ_STATISTIC\t%.5lg\n",2*(est.llk1-est.llk0));
  if ( fixAlpha < 0 )
    hprintf(wf, "SAME_ANCESTRY\tESTIMATED_CONTAM\t%.5lg\n",alpha);
  if ( fixPC.empty() ) {
    for(int32_t i=0; i < numPC; ++i)   
      hprintf(wf, "SAME_ANCESTRY\tINTENDED_PC%d\t%.5lg\n",i+1,pc1[i]);
  }

  bool use_diff_model = false;
  
  if ( ( fixPC.empty() ) && ( !withinAncestry ) ) {
    double alpha_b = alpha;
    double* pc1_b = new double[numPC];
    double* pc2_b = new double[numPC];
    
    for(int32_t i=0; i < numPC; ++i) {
      pc1_b[i] = pc1[i];
      pc2_b[i] = pc2[i];
    }
    est.optimizeLLK(&alpha_b, pc1_b, pc2_b, false);
  
    hprintf(wf, "DIFF_ANCESTRY\tLOG_LIKELIHOOD_H0\t%.5lg\n",est.llk0);
    hprintf(wf, "DIFF_ANCESTRY\tLOG_LIKELIHOOD_H1\t%.5lg\n",est.llk1);
    hprintf(wf, "DIFF_ANCESTRY\tCHISQ_STATISTIC\t%.5lg\n",2*(est.llk1-est.llk0));
    if ( fixAlpha < 0 )
      hprintf(wf, "DIFF_ANCESTRY\tESTIMATED_CONTAM\t%.5lg\n",alpha_b);
    if ( fixPC.empty() ) {
      for(int32_t i=0; i < numPC; ++i)   
	hprintf(wf, "DIFF_ANCESTRY\tINTENDED_PC%d\t%.5lg\n",i+1,pc1_b[i]);
      for(int32_t i=0; i < numPC; ++i)   
	hprintf(wf, "DIFF_ANCESTRY\tCONTAM_PC%d\t%.5lg\n",i+1,pc2_b[i]);      
    }

    if ( 2*(est.llk1-est.llk0) > 2.70 ) {
      alpha = alpha_b;
      for(int32_t i=0; i < numPC; ++i) {
	pc1[i] = pc1_b[i];
	pc2[i] = pc2_b[i];
      }
      use_diff_model = true;
      llk1 = est.llk1;
    }
    delete [] pc1_b;
    delete [] pc2_b;
  }

  double probRef[3] = {1.0, 0.5, 0.0};

  if ( refBias ) {
    double alpha_b = alpha;
    double* pc1_b = new double[numPC];
    double* pc2_b = new double[numPC];

    for(int32_t i=0; i < numPC; ++i) {
      pc1_b[i] = pc1[i];
      pc2_b[i] = pc2[i];
    }
    est.optimizeLLK(&alpha_b, pc1_b, pc2_b, withinAncestry, true);
  
    hprintf(wf, "REF_BIAS\tLOG_LIKELIHOOD_H0\t%.5lg\n",est.llk0);
    hprintf(wf, "REF_BIAS\tLOG_LIKELIHOOD_H1\t%.5lg\n",est.llk1);
    hprintf(wf, "REF_BIAS\tCHISQ_STATISTIC\t%.5lg\n",2*(est.llk1-est.llk0));
    if ( fixAlpha < 0 )
      hprintf(wf, "REF_BIAS\tESTIMATED_CONTAM\t%.5lg\n",alpha_b);
    if ( fixPC.empty() ) {
      for(int32_t i=0; i < numPC; ++i)   
	hprintf(wf, "REF_BIAS\tINTENDED_PC%d\t%.5lg\n",i+1,pc1_b[i]);
      for(int32_t i=0; i < numPC; ++i)   
	hprintf(wf, "REF_BIAS\tCONTAM_PC%d\t%.5lg\n",i+1,pc2_b[i]);      
    }
    hprintf(wf, "REF_BIAS\tPROB_REF_RR\t%.5lg\n",est.probRef[0]);
    hprintf(wf, "REF_BIAS\tPROB_REF_RA\t%.5lg\n",est.probRef[1]);
    hprintf(wf, "REF_BIAS\tPROB_REF_AA\t%.5lg\n",est.probRef[2]);

    if ( 2*(est.llk1-est.llk0) > 2.70 ) {
      alpha = alpha_b;
      for(int32_t i=0; i < numPC; ++i) {
	pc1[i] = pc1_b[i];
	pc2[i] = pc2_b[i];
      }
      use_diff_model = true;
      llk1 = est.llk1;
      
      probRef[0] = est.probRef[0];
      probRef[1] = est.probRef[1];
      probRef[2] = est.probRef[2];      
    }
    delete [] pc1_b;
    delete [] pc2_b;    
  }

  hprintf(wf, "FINAL_MODEL\tUSE_DIFF_MODEL\t%d\n",use_diff_model ? 1 : 0);  
  hprintf(wf, "FINAL_MODEL\tLOG_LIKELIHOOD_H0\t%.5lg\n",llk0);
  hprintf(wf, "FINAL_MODEL\tLOG_LIKELIHOOD_H1\t%.5lg\n",llk1);
  hprintf(wf, "FINAL_MODEL\tCHISQ_STATISTIC\t%.5lg\n",2*(llk1-llk0));
  if ( fixAlpha < 0 )
      hprintf(wf, "FINAL_MODEL\tESTIMATED_CONTAM\t%.5lg\n",alpha);
  if ( fixPC.empty() ) {
    for(int32_t i=0; i < numPC; ++i)   
      hprintf(wf, "FINAL_MODEL\tINTENDED_PC%d\t%.5lg\n",i+1,pc1[i]);
    for(int32_t i=0; i < numPC; ++i)   
      hprintf(wf, "FINAL_MODEL\tCONTAM_PC%d\t%.5lg\n",i+1,pc2[i]);      
  }
  hprintf(wf, "FINAL_MODEL\tPROB_REF_RR\t%.5lg\n",probRef[0]);
  hprintf(wf, "FINAL_MODEL\tPROB_REF_RA\t%.5lg\n",probRef[1]);
  hprintf(wf, "FINAL_MODEL\tPROB_REF_AA\t%.5lg\n",probRef[2]);  

  /*
  double* rdps = est.get_relative_depths(pc1);
  hprintf(wf, "RELATIVE_DEPTH\tOR_REF\t%.5lg\n", rdps[0]);
  hprintf(wf, "RELATIVE_DEPTH\tOR_HET\t%.5lg\n", rdps[1]);
  hprintf(wf, "RELATIVE_DEPTH\tOR_ALT\t%.5lg\n", rdps[2]);
  hprintf(wf, "RELATIVE_DEPTH\tOBS_REF_REF\t%.5lg\n", rdps[3]);
  hprintf(wf, "RELATIVE_DEPTH\tOBS_REF_ALT\t%.5lg\n", rdps[4]);
  hprintf(wf, "RELATIVE_DEPTH\tOBS_HET_REF\t%.5lg\n", rdps[5]);
  hprintf(wf, "RELATIVE_DEPTH\tOBS_HET_ALT\t%.5lg\n", rdps[6]);
  hprintf(wf, "RELATIVE_DEPTH\tOBS_ALT_REF\t%.5lg\n", rdps[7]);
  hprintf(wf, "RELATIVE_DEPTH\tOBS_ALT_ALT\t%.5lg\n", rdps[8]);
  hprintf(wf, "RELATIVE_DEPTH\tEXP_REF_REF\t%.5lg\n", rdps[9]);
  hprintf(wf, "RELATIVE_DEPTH\tEXP_REF_ALT\t%.5lg\n", rdps[10]);
  hprintf(wf, "RELATIVE_DEPTH\tEXP_HET_REF\t%.5lg\n", rdps[11]);
  hprintf(wf, "RELATIVE_DEPTH\tEXP_HET_ALT\t%.5lg\n", rdps[12]);
  hprintf(wf, "RELATIVE_DEPTH\tEXP_ALT_REF\t%.5lg\n", rdps[13]);
  hprintf(wf, "RELATIVE_DEPTH\tEXP_ALT_ALT\t%.5lg\n", rdps[14]);  
  delete [] rdps;
  */
  
  hts_close(wf);
  
  if ( writePileup ) {
    htsFile* plpf = hts_open( (outPrefix + ".plp.gz").c_str(), "wz" );
    if ( plpf == NULL )
      error("Cannot open %s for writing", (outPrefix+".plp.gz").c_str());
    est.writePileup(plpf, alpha, pc1, pc2);
    hts_close(plpf);
  }
  

  delete [] pc1;
  delete [] pc2;
  delete [] est.fixPC;

  return 0;
}

