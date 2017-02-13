#include "cramore.h"
#include "bcf_ordered_reader.h"

int32_t cmdCramVerifyPairID(int32_t argc, char** argv) {
  std::string inSam; // SAM, BAM, or CRAM
  std::string inVcf; // VCF or VCF with GT fields
  std::string outPrefix;
  std::string tagGroup;
  std::string tagUMI;
  int32_t capBQ = 40;
  int32_t minBQ = 13;
  int32_t minMQ = 20;
  int32_t minTD = 0;
  int32_t qcExclFlag = 0x0f04;
  std::vector<std::string> smIDs;
  std::vector<double> gridAlpha;
  std::vector<double> gridASE;
  int32_t verbose = 10000;
  //int32_t minMAC = 1;
  //int32_t minCallRate = 0.9;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Options for input SAM/BAM/CRAM", NULL)
    LONG_STRING_PARAM("sam",&inSam, "Input SAM/BAM/CRAM file. Must be sorted by coordinates and indexed")
    LONG_STRING_PARAM("tag-group",&tagGroup, "Tag representing readgroup or cell barcodes, in the case to partition the BAM file into multiple groups. For 10x genomics, use CB")
    LONG_STRING_PARAM("tag-UMI",&tagUMI, "Tag representing UMIs. For 10x genomiucs, use UB")

    LONG_PARAM_GROUP("Options for input VCF/BCF", NULL)
    LONG_STRING_PARAM("vcf",&inVcf, "Input VCF/BCF file, containing the individual genotypes as GT field")
    LONG_MULTI_STRING_PARAM("sm",&smIDs, "List of sample IDs to compare to (default: use all)")    

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out",&outPrefix,"Output file prefix")

    LONG_PARAM_GROUP("Analysis Options", NULL)
    LONG_INT_PARAM("cap-BQ", &capBQ, "Maximum base quality (higher BQ will be capped)")
    LONG_INT_PARAM("min-BQ", &minBQ, "Minimum base quality to consider (lower BQ will be skipped)")
    LONG_INT_PARAM("min-MQ", &minMQ, "Minimum mapping quality to consider (lower MQ will be ignored)")
    LONG_INT_PARAM("min-TD", &minTD, "Minimum distance to the tail (lower will be ignored)")
    LONG_INT_PARAM("excl-flag", &qcExclFlag, "SAM/BAM FLAGs to be excluded")    
    LONG_MULTI_DOUBLE_PARAM("alpha",&gridAlpha, "Grid of alpha to search for (default is 0, 0.1, 0.2, 0.3, 0.4  0.5)")
    LONG_MULTI_INT_PARAM("ase",&gridASE, "Grid of allele-specific expression to search for (default is 0, 0.2, 0.3, 0.4, 0.5) -- Not implemented")
    LONG_INT_PARAM("verbose",&verbose, "Verbose message frequency")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  if ( gridAlpha.empty() ) {
    gridAlpha.push_back(0);    
    gridAlpha.push_back(0.02);
    gridAlpha.push_back(0.05);
    gridAlpha.push_back(0.1);
    gridAlpha.push_back(0.15);
    gridAlpha.push_back(0.2);
    gridAlpha.push_back(0.25);        
    gridAlpha.push_back(0.3);
    gridAlpha.push_back(0.35);    
    gridAlpha.push_back(0.4);
    gridAlpha.push_back(0.45);    
    gridAlpha.push_back(0.5);    
  }

  if ( gridASE.empty() ) {
    gridASE.push_back(0.1);
    gridASE.push_back(0.2);
    gridASE.push_back(0.3);
    gridASE.push_back(0.4);
    gridASE.push_back(0.5);    
  }

  // load VCF files. This VCF should only contain hard genotypes in GT field
  std::vector<GenomeInterval> intervals;    
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();

  samFile* in = NULL;
  bam_hdr_t *header = NULL;

  if ( ( in = sam_open(inSam.c_str(), "r") ) == 0 ) {
    error("[E:%s:%d %s] Cannot open file %s\n",__FILE__,__LINE__,__FUNCTION__,inSam.c_str());    
  }

  if ( ( header = sam_hdr_read(in) ) == 0 ) {
    error("[E:%s:%d %s] Cannot open header from %s\n",__FILE__,__LINE__,__FUNCTION__,inSam.c_str());
  }

  if ( outPrefix.empty() )
    error("[E:%s:%d %s] --out parameter is missing",__FILE__,__LINE__,__FUNCTION__);

  bam1_t *b = bam_init1();
  //int32_t r;    

  int32_t nsamples = bcf_hdr_nsamples(odr.hdr);
  std::vector< std::vector<uint8_t> > v_gts;
  std::vector<double> v_afs;
  std::vector<int32_t> v_rids;
  std::vector<int32_t> v_poss;
  std::vector<char> v_refs;
  std::vector<char> v_alts;
  
  // identify samples to focus on
  std::vector<int> vids;
  std::vector<std::string> vSMs;
  if ( smIDs.empty() ) {
    for(int32_t i=0; i < nsamples; ++i) {
      vids.push_back(i);
      vSMs.push_back(odr.hdr->samples[i]);
    }
  }
  else {
    std::set<std::string> smMap;
    for(int32_t i=0; i < (int32_t)smIDs.size(); ++i) {
      smMap.insert(smIDs[i]);
    }
    for(int32_t i=0; i < nsamples; ++i) {
      if ( smMap.find(odr.hdr->samples[i]) != smMap.end() ) {
	vids.push_back(i);
	vSMs.push_back(odr.hdr->samples[i]);
      }
    }
  }
  int32_t nv = (int32_t)vids.size();

  notice("Started reading from the VCF file %s", inVcf.c_str());
  
  // read VCF and store genotypes
  while( odr.read(iv) ) { // read marker
    //bool skip = false;
    bcf_unpack(iv, BCF_UN_ALL);
    
    if ( iv->n_allele > 2 ) continue; // skip multi-allelics
    if ( !bcf_is_snp(iv) ) continue;  // focus only on SNPs
    
    // chrom is iv->rid
    // position is iv->pos
    // ref is iv->d.allele[0]
    // read genotypes
    int32_t rid = iv->rid;
    int32_t pos = iv->pos;
    char ref = iv->d.allele[0][0];
    char alt = iv->d.allele[1][0];

    // store marker information
    v_rids.push_back(rid);
    v_poss.push_back(pos);
    v_refs.push_back(ref);
    v_alts.push_back(alt);
    
    uint32_t* p_gts = (uint32_t*)calloc(nsamples * 2, sizeof(uint32_t));
    int32_t n_gts = 0;
    
    // extract genotypes fpr selected individuals
    if ( bcf_get_genotypes(odr.hdr, iv, &p_gts, &n_gts) < 0 ) {
      error("[E:%s:%d %s] Cannot extract genotypes at %s:%d %c/%c",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, rid), pos, ref, alt);
    }
    if ( n_gts != nsamples * 2 ) {
      error("[E:%s:%d %s] Cannot extract %d genotypes at %s:%d %c/%c. Extracted only %d",__FILE__,__LINE__,__FUNCTION__, nsamples * 2, bcf_hdr_id2name(odr.hdr, rid), pos, ref, alt, n_gts); 
    }

    v_gts.resize( v_gts.size() + 1 );
    //v_gts.push_back( std::vector<uint8_t>(nv, 0) );
    std::vector<uint8_t>& v = v_gts.back();
    v.resize(nv);

    int32_t ac = 0, an = 0;
    for(int32_t i=0; i < nv; ++i) {   // bi-allelic encoding of variant
      int32_t g1 = p_gts[2*vids[i]];
      int32_t g2 = p_gts[2*vids[i]+1];
      uint8_t geno;
      if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
	geno = 0;
      }
      else {
	geno = (uint8_t)(((bcf_gt_allele(g1) > 0) ? 1 : 0) + ((bcf_gt_allele(g2) > 0) ? 1 : 0) + 1);
	ac += (geno-1);
	an += 2;
      }
      v[i] = geno;
    }

    if ( ( verbose > 0 ) && ( v_rids.size() % verbose == 0 ) )
      notice("Reading %d markers from the VCF file", (int32_t)v_rids.size());

    v_afs.push_back((double)(ac+5e-7)/(double)(an+1e-6));

    free(p_gts);
  }

  notice("Finished reading %d markers from the VCF file", (int32_t)v_rids.size());  

  // calculate genotype likelihoods from BAM/CRAMs
  // we expect
  // Pr(B|AA,AA) = (e/3) ~ e
  // Pr(B|AA,AB) = a(1/2-e/6) ~ a/2
  // Pr(B|AA,BB) = a(1-e) ~ a
  // Pr(B|AB,AA) = (1-a)(1/2-e/6) ~ (1-a)/2
  // Pr(B|AB,AB) = (1/2-e/6) ~ 1/2
  // Pr(B|AB,BB) = (1-a)(1/2-e/6)+(1-e)a ~ (1-a)/2+a = (1+a)/2
  // Pr(B|BB,AA) = (1-a)(1-e) ~ (1-a)
  // Pr(B|BB,AB) = (1-a)(1-e)+a(1/2-e/t) ~ 1-a/2
  // Pr(B|BB,BB) = (1-e) ~ 1
  // llk for [n * n * m * a * bc] ??  0.1 0.3 0.5 0.7 0.9
  
  int32_t nAlpha = (int32_t)gridAlpha.size();

  kstring_t readseq = {0,0,0};
  kstring_t readqual = {0,0,0};	
  
  hts_idx_t *idx = sam_index_load(in, inSam.c_str());
  if ( idx == NULL )
    error("[E:%s:%d %s] Cannot load index file for %s",__FILE__,__LINE__,__FUNCTION__,inSam.c_str());

  std::vector< std::vector<uint8_t> > bqs;
  std::vector< std::vector<uint8_t> > alleles; // 0 - REF, 1 - ALT, 2 - OTHER
  std::vector< std::vector<uint32_t> > ibcds;
  std::map< std::string, uint32_t > bcMap;
  std::vector< int32_t > bcReads(1,0);
  std::vector< int32_t > bcVars(1,0);  
  bcMap["."] = 0;
  
  char reg[255];
  char gtag[2] = {0,0};
  char utag[2] = {0,0};  

  int32_t nReadsAll = 0, nReadsPass = 0, nReadsRedundant = 0, nReadsN = 0, nReadsLQ = 0;

  if ( tagGroup.empty() ) { // do nothing
  }
  else if ( tagGroup.size() == 2 ) {
    gtag[0] = tagGroup.at(0);
    gtag[1] = tagGroup.at(1);    
  }
  else {
    error("[E:%s:%d %s] Cannot recognize group tag %s. It is suppose to be a length 2 string",__FILE__,__LINE__,__FUNCTION__,tagGroup.c_str());
  }

  if ( tagUMI.empty() ) { // do nothing
  }
  else if ( tagUMI.size() == 2 ) {
    utag[0] = tagUMI.at(0);
    utag[1] = tagUMI.at(1);    
  }
  else {
    error("[E:%s:%d %s] Cannot recognize UMI tag %s. It is suppose to be a length 2 string",__FILE__,__LINE__,__FUNCTION__,tagUMI.c_str());
  }  
  
  for(int32_t j=0; j < (int32_t)v_poss.size(); ++j) {  // process each marker separately
    if ( (j+1) % 10000 == 0 )
      notice("Finished Processing %d reads across %d variants across %d barcodes, filtering %d (%.2lf%%) reads, including %d (%.2lf%%) gapped reads, %d (%.2lf%%) low quality reads, and %d (%.2lf%%) redundant/qcfail reads from the BAM file %s", nReadsPass, j+1, (int32_t)bcMap.size(), nReadsAll - nReadsPass, 100.0 * (nReadsAll - nReadsPass) / nReadsAll, nReadsN, 100.0 * nReadsN / nReadsAll, nReadsLQ, 100.0 * nReadsLQ / nReadsAll, nReadsRedundant, 100.0 * nReadsRedundant / nReadsAll,  inSam.c_str());            
    
    bqs.resize( bqs.size() + 1 );
    alleles.resize( alleles.size() + 1 );
    ibcds.resize( ibcds.size() + 1 );

    sprintf(reg, "%s:%d-%d", bcf_hdr_id2name(odr.hdr, v_rids[j]), v_poss[j], v_poss[j]);

    std::vector<uint8_t>& v_bq = bqs.back();
    std::vector<uint8_t>& v_al = alleles.back();
    std::vector<uint32_t>& v_ibcd = ibcds.back();
    std::set<std::string> sUMI;    
    
    char ref = v_refs[j];
    char alt = v_alts[j];
    char base, qual;
    int32_t rpos;

    std::set<int32_t> s_ibcds;

    //std::string sal;
    //sal += ref;
    //sal += alt;
    //sal += "/";
    
    hts_itr_t* itr = bam_itr_querys(idx, header, reg);
    while( sam_itr_next(in, itr, b) >= 0 ) {
      ++nReadsAll;

      //fprintf(stderr,"\n%d ",b->core.qual);

      if ( b->core.flag & qcExclFlag ) {
	++nReadsRedundant;	
	continue;
      }
      
      if ( b->core.qual < minMQ ) {
	++nReadsN;
	continue;
      }
           
 
      bam_get_base_and_qual_and_read_and_qual(b, (uint32_t)v_poss[j], base, qual, rpos, &readseq, &readqual);
      
      if ( base == 'N' ) {
	++nReadsN;
	continue;
      }

      if ( qual-33 < minBQ ) { ++nReadsLQ; continue; }
      if ( rpos < minTD-1 ) { ++nReadsLQ; continue; }
      if ( rpos + minTD > b->core.l_qseq ) { ++nReadsLQ; continue; }

      uint8_t *bcd = (*gtag) ? (uint8_t*) bam_aux_get(b, gtag) : NULL;
      uint32_t ibcd = 0;
      const char* sbcd = ".";
      if ( bcd != NULL ) {
	if ( *bcd == 'Z' ) {
	  sbcd = bam_aux2Z(bcd);
	  if ( bcMap.find(sbcd) == bcMap.end() ) {
	    ibcd = bcMap.size();
	    bcMap[sbcd] = ibcd;
	    bcReads.resize(bcMap.size(),0);
	    bcVars.resize(bcMap.size(),0);	    
	  }
	  else {
	    ibcd = bcMap[sbcd];
	  }
	}
      }

      uint8_t *umi = (*utag) ? (uint8_t*) bam_aux_get(b, utag) : NULL;
      if ( umi != NULL ) {
	if ( *umi == 'Z' ) {
	  char* sumi = bam_aux2Z(umi);
	  std::string bcumi = std::string(sbcd)+sumi;
	  if ( sUMI.find(sumi) == sUMI.end() ) {
	    sUMI.insert(sumi);
	  }
	  else {
	    ++nReadsRedundant;
	    continue; // skip UMI already seen at the same position
	  }
	}
      }

      v_al.push_back( ( base == ref ) ? 0 : ( ( base == alt ) ? 1 : 2 ) );
      v_bq.push_back( qual-33 > capBQ ? capBQ : qual-33 );
      v_ibcd.push_back( ibcd );
      ++bcReads[ibcd];
      s_ibcds.insert(ibcd);

      //sal += (( base == ref ) ? "0" : ( ( base == alt ) ? "1" : "2" ) );
    //sal += qual;
    //sal += (char)(33+rpos);
    //sal += ",";
      
      ++nReadsPass;
    }
    sam_itr_destroy(itr);

    //if ( !v_al.empty() ) 
    //notice("%s",sal.c_str());

    for(std::set<int32_t>::const_iterator it = s_ibcds.begin(); it != s_ibcds.end(); ++it) {
      ++bcVars[*it];
    }
  }

  odr.close();

  notice("Finished processing %d reads across %d variants across %d barcodes, filtering %d (%.2lf%%) reads, including %d (%.2lf%%) gapped reads, %d (%.2lf%%) low quality reads, and %d (%.2lf%%) redundant/qcfail reads from the BAM file %s", nReadsPass, (int32_t)v_poss.size(), (int32_t)bcMap.size(), nReadsAll - nReadsPass, 100.0 * (nReadsAll - nReadsPass) / nReadsAll, nReadsN, 100.0 * nReadsN / nReadsAll, nReadsLQ, 100.0 * nReadsLQ / nReadsAll, nReadsRedundant, 100.0 * nReadsRedundant / nReadsAll,  inSam.c_str());      

  // start evaluating genotype concordances
  // calculate for (nBcd) x (nInds) to find the best matching genotypes first
  notice("Starting to identify best matching individual IDs");

  htsFile* wsingle = hts_open((outPrefix+".single").c_str(),"w");
  htsFile* wpair = hts_open((outPrefix+".pair").c_str(),"w");
  htsFile* wbest = hts_open((outPrefix+".best").c_str(),"w");

  hprintf(wsingle, "BARCODE\tSM_ID\tN_READS\tN_SNPS\tLLK1\tLLK0\n");
  hprintf(wpair,   "BARCODE\tSM1_ID\tSM2_ID\tALPHA\tN_READS\tN_SNPS\tLLK12\tLLK1\tLLK0\tLLK10\tLLK00\n");
  hprintf(wbest,   "BARCODE\tSM1_ID\tSM2_ID\tALPHA\tN_READS\tN_SNPS\tLLK12\tLLK1\tLLK0\tLLK10\tLLK00\n");    

  if ( ( wbest == NULL ) || ( wpair == NULL ) || ( wsingle == NULL ) )
    error("[E:%s:%d %s] Cannot create %s.single, %s.pair, and %s.best files",__FILE__,__LINE__,__FUNCTION__,outPrefix.c_str(), outPrefix.c_str());
  
  int32_t nbcd = (int32_t)bcMap.size();
  std::vector<double> llks(nbcd * nv, 0);
  std::vector<double> llk0s(nbcd, 0);  
  int32_t nsnp = (int32_t)bqs.size();
  for(int32_t i=0; i < nsnp; ++i) {
    if ( ( verbose > 0 ) && ( (i+1) % verbose == 0 ) )
      notice("Processing %d markers...",i+1);

    if ( bqs[i].size() == 0 ) continue;
    
    std::vector<double> GLs(nbcd * 4, 0);

    std::vector<uint8_t>& vgt = v_gts[i];
    
    for(int32_t j=0; j < (int32_t)bqs[i].size(); ++j) {
      uint8_t bq = bqs[i][j];
      uint8_t al = alleles[i][j];
      uint32_t ibcd = ibcds[i][j];

      if ( al == 2 ) continue;

      GLs[ibcd * 4 + 1] += ( (al == 0) ? phredConv.phred2LogMat3[bq] : -0.1*bq );
      GLs[ibcd * 4 + 2] += phredConv.phred2HalfLogMat3[bq];
      GLs[ibcd * 4 + 3] += ( (al == 1) ? phredConv.phred2LogMat3[bq] : -0.1*bq );      
    }

    double p0 = (1.-v_afs[i])*(1.-v_afs[i]);
    double p1 = 2. * v_afs[i] * (1.-v_afs[i]);
    double p2 = v_afs[i] * v_afs[i];

    // calculate genotype likelihoods for each batcodes
    for(int32_t j=0; j < nbcd; ++j) {
      GLs[j * 4] = ( GLs[j * 4 + 1] > GLs[j * 4 + 2] ) ? ( GLs[j * 4 + 1] > GLs[j * 4 + 3] ? GLs[j * 4 + 1] : GLs[j * 4 + 3] ) : ( GLs[j * 4 + 2] > GLs[j * 4 + 3] ? GLs[j * 4 + 2] : GLs[j * 4 + 3] );
      GLs[j*4+1] -= GLs[j*4];
      GLs[j*4+2] -= GLs[j*4];
      GLs[j*4+3] -= GLs[j*4];
      GLs[j*4] = log10(p0*pow(10.,GLs[j*4+1]) + p1*pow(10.,GLs[j*4+2]) + p2*pow(10.,GLs[j*4+3]));

      for(int32_t k=0; k < nv; ++k) {
	llks[j * nv + k] += GLs[j*4 + vgt[k]];
      }
      llk0s[j] += GLs[j*4];
    }
  }

  std::vector<std::string> sbcd(bcMap.size());
  for(std::map<std::string,uint32_t>::const_iterator it = bcMap.begin(); it != bcMap.end(); ++it) {
    sbcd[it->second] = it->first;
  }

  // find the best matching individual
  std::vector<int32_t> iBest(nbcd,0);
  std::vector<double> llkBest(nbcd,0);
  notice("Identifying best-matching individual..");
  for(int32_t i=0; i < nbcd; ++i) {
    double imax = -1;
    double maxLLK = -1e300;
    //double inext = -1;
    double nextLLK = -1e300;
    for(int32_t j=0; j < nv; ++j) {
      hprintf(wsingle,"%s\t%s\t%d\t%d\t%.5lf\t%.5lf\n", sbcd[i].c_str(), vSMs[j].c_str(), bcReads[i], bcVars[i], llks[i*nv+j], llk0s[i]);
      
      if ( llks[i*nv + j] > maxLLK ) {
	//inext = imax;
	nextLLK = maxLLK;
	imax = j;
	maxLLK = llks[i*nv+j];
      }
      else if ( llks[i*nv + j] > nextLLK ) {
	//inext = j;
	nextLLK = llks[i*nv+j];
      }
    }

    //notice("Barcode:%s\tnReads:%d\tnVariants:%d\tBest:%s\tNext:%s\tmaxLLK=%.5lf\tnextLLK=%.5lf\t%LLK0=%.5lf", sbcd[i].c_str(), bcReads[i], bcVars[i], vSMs[imax].c_str(), vSMs[inext].c_str(), llks[i*nv+imax], nextLLK, llk0s[i]);
    iBest[i] = imax;
    llkBest[i] = maxLLK;
  }

  hts_close(wsingle);

  // start finding the next-best matching individual
  std::vector<double> pllks(nAlpha * nbcd * nv, 0);
  std::vector<double> pllk0s(nAlpha * nbcd , 0);
  std::vector<double> pllk00s(nAlpha * nbcd , 0);  
  
  for(int32_t i=0; i < nsnp; ++i) {
    if ( ( verbose > 0 ) && ( (i+1) % verbose == 0 ) )
      notice("Processing %d markers...",i+1);

    if ( bqs[i].size() == 0 ) continue;
    
    std::vector<double> pGs(nbcd * nAlpha * 16, 1.);
    std::vector<uint8_t>& vgt = v_gts[i];
    std::vector<int32_t> cbcds(nbcd,0);
    for(int32_t j=0; j < (int32_t)bqs[i].size(); ++j) {
      uint8_t bq = bqs[i][j];
      uint8_t al = alleles[i][j];
      uint32_t ibcd = ibcds[i][j];

      if ( al == 2 ) continue;

      double pR = (al == 0) ? phredConv.phred2Mat3[bq] : phredConv.phred2Err[bq];
      double pA = (al == 1) ? phredConv.phred2Mat3[bq] : phredConv.phred2Err[bq];

      double maxpG = 0;      
      for(int32_t k=0; k < nAlpha; ++k) {
	for(int32_t l=0; l < 3; ++l) {  // 1-Alpha
	  for(int32_t m=0; m < 3; ++m) { // Alpha
	    double p = 0.5*l + (m-l)*0.5*gridAlpha[k]; // %A
	    double& pG = pGs[ibcd*16*nAlpha + k*16 + (l+1)*4 + m+1];
	    pG *= (pR * (1-p) + pA * p);
	    if ( maxpG < pG )
	      maxpG = pG;
	  }
	}
      }

      for(int32_t k=0; k < nAlpha; ++k) {
	// normalize
	for(int32_t l=0; l < 3; ++l) {  // 1-Alpha
	  for(int32_t m=0; m < 3; ++m) { // Alpha
	    pGs[ibcd*16*nAlpha + k*16 + (l+1)*4 + m+1] /= maxpG;
	  }
	}
      }
      ++cbcds[ibcd];
    }

    double p0 = (1.-v_afs[i])*(1.-v_afs[i]);
    double p1 = 2.*v_afs[i]*(1.-v_afs[i]);
    double p2 = v_afs[i] * v_afs[i];
    for(int32_t j=0; j < nbcd; ++j) {
      if ( cbcds[j] == 0 ) continue;
      for(int32_t k=0; k < nAlpha; ++k) {
	for(int32_t l=0; l < 3; ++l) {
	  pGs[j*16*nAlpha + k*16 + (l+1)*4] = p0*pGs[j*16*nAlpha + k*16 + (l+1)*4 +1] + p1*pGs[j*16*nAlpha + k*16 + (l+1)*4 +2] + p2*pGs[j*16*nAlpha + k*16 + (l+1)*4 +2];
	  pGs[j*16*nAlpha + k*16 + (l+1)] = p0*pGs[j*16*nAlpha + k*16 + l + 5] + p1*pGs[j*16*nAlpha + k*16 + l + 9] + p2*pGs[j*16*nAlpha + k*16 + l + 13];
	}
	pGs[j*16*nAlpha + k*16] = p0*p0*pGs[j*16*nAlpha + k*16 + 5] + p0*p1*pGs[j*16*nAlpha + k*16 + 6] + p0*p2*pGs[j*16*nAlpha + k*16 + 7] + p1*p0*pGs[j*16*nAlpha + k*16 + 9] + p1*p1*pGs[j*16*nAlpha + k*16 + 10] + p1*p2*pGs[j*16*nAlpha + k*16 + 11] + p2*p0*pGs[j*16*nAlpha + k*16 + 13] + p2*p1*pGs[j*16*nAlpha + k*16 + 14] + p2*p2*pGs[j*16*nAlpha + k*16 + 15];
	
	for(int32_t l=0; l < 16; ++l) {
	  pGs[j*16*nAlpha + k*16 + l] = log10(pGs[j*16*nAlpha + k*16 + l]);
	}

	for(int32_t l=0; l < nv; ++l) {
	  pllks[j*nAlpha*nv + k*nv + l] += pGs[j*16*nAlpha + k*16 + vgt[iBest[j]]*4 + vgt[l]];
	}
	pllk0s[j*nAlpha + k] += pGs[j*16*nAlpha + k*16 + vgt[iBest[j]]*4];
	pllk00s[j*nAlpha + k] += pGs[j*16*nAlpha + k*16];	
      }
    }
  }

  for(int32_t i=0; i < nbcd; ++i) {
    double jBest = 0;
    double kBest = 0;
    double maxLLK = -1e300;
    for(int32_t j=0; j < nAlpha; ++j) {
      for(int32_t k=0; k < nv; ++k) {
	hprintf(wpair,"%s\t%s\t%s\t%.3lf\t%d\t%d\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\n", sbcd[i].c_str(), vSMs[iBest[i]].c_str(), vSMs[k].c_str(), gridAlpha[j], bcReads[i], bcVars[i], pllks[i*nAlpha*nv + j*nv +k], pllks[i*nAlpha*nv], pllk00s[i*nAlpha], pllk0s[i*nAlpha + j], pllk00s[i*nAlpha + j]);
	
	if ( pllks[i*nAlpha*nv + j*nv + k] > maxLLK ) {
	  jBest = j;
	  kBest = k;
	  maxLLK = pllks[i*nAlpha*nv + j*nv + k];
	}
      }
    }

    //notice("Barcode: %s\tBest: (%s, %s, %.1lf%%)\tmaxLLK=%.5lf", sbcd[i].c_str(), vSMs[iBest[i]].c_str(), vSMs[kBest].c_str(), 100*gridAlpha[jBest], maxLLK);
    hprintf(wbest,"%s\t%s\t%s\t%.3lf\t%d\t%d\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\n", sbcd[i].c_str(), vSMs[iBest[i]].c_str(), vSMs[kBest].c_str(), gridAlpha[jBest], bcReads[i], bcVars[i], pllks[i*nAlpha*nv + jBest*nv +kBest], pllks[i*nAlpha*nv], pllk00s[i*nAlpha], pllk0s[i*nAlpha + jBest], pllk00s[i*nAlpha + jBest]);    
  }

  notice("Finished writing output files");

  hts_close(wpair);
  hts_close(wbest);
  
  return 0;
}
