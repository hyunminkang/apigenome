#include "cramore.h"
#include "estimator.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"

int32_t cmdCramSparseGenotype(int32_t argc, char** argv) {
  std::string inVcf;
  std::vector<std::string> inCrams;
  std::string inCramList;
  std::string out;
  std::string reg;
  int32_t threads = 1;
  double gl_adj = 0.01;
  int32_t xStart = 2699520;
  int32_t xStop = 154931044;
  std::string xLabel("X");
  std::string yLabel("Y");
  std::string mtLabel("MT");
  std::string sexMap;

  #ifdef _OPENMP
  threads = 4;
  #endif
  
  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

    LONG_PARAM_GROUP("Input Sequence Files", NULL)
    LONG_MULTI_STRING_PARAM("in-cram",&inCrams, "Input CRAM file(s)")
    LONG_STRING_PARAM("in-cram-list",&inCramList, "File containing input CRAM files")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file name")

    LONG_PARAM_GROUP("Sex Chromosomes",NULL)
    LONG_STRING_PARAM("xLabel", &xLabel, "Contig name for X chromosome")
    LONG_STRING_PARAM("yLabel", &xLabel, "Contig name for Y chromosome")
    LONG_STRING_PARAM("mtLabel", &xLabel, "Contig name for MT chromosome")
    LONG_INT_PARAM("xStart", &xStart, "Start base position of non-PAR region in X chromosome")
    LONG_INT_PARAM("xStop",  &xStop,  "End base position of non-PAR region in X chromosome")    

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("threads",&threads, "Number of threads to parallelize")
    LONG_STRING_PARAM("sex-map",&sexMap, "Sex map file, containing ID and sex (1 for male and 2 for female) for each individual")
    LONG_DOUBLE_PARAM("gl-adj",&gl_adj, "Genotype likelihood adjustment factor at homozygous sites as Pr(1|0/0) or Pr(0|1/1)")
  END_LONG_PARAMS();
  
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();
  
  // sanity check of input arguments
  if ( inVcf.empty() || out.empty() || ( inCrams.empty() && inCramList.empty() ) ) {
    error("[E:%s:%d %s] --in-vcf, --out, --in-cram (or --in-cram-list) are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  if ( ( !inCrams.empty() ) && ( ! inCramList.empty() ) ) {
    error("[E:%s:%d %s] --in-cram-list and --in-cram options cannot be used together",__FILE__,__LINE__,__FUNCTION__);
  }

  if ( ! inCramList.empty() ) {
    htsFile* fp = hts_open(inCramList.c_str(),"r");
    if ( fp == NULL )
      error("[E:%s:%d %s] Cannot open file %s for reading",__FILE__,__LINE__,__FUNCTION__, inCramList.c_str());


    int32_t lstr = 0;
    int32_t* fields = NULL;
    int32_t n = 0;
    kstring_t str = {0,0,0};      

    while( ( lstr = hts_getline(fp, KS_SEP_LINE, &str) ) >= 0 ) {
      fields = ksplit(&str, 0, &n);
      if ( n > 1 )
	error("[E:%s:%d %s] in-cram-list file %s contains whitespace - # fields = %d, (%s, %s)",__FILE__,__LINE__,__FUNCTION__, inCramList.c_str(), n, str.s + fields[0], str.s + fields[1]);
      inCrams.push_back(std::string(str.s));
    }
    hts_close(fp);
  }

  std::vector<GenomeInterval> intervals;
  if ( !reg.empty() ) {
    parse_intervals(intervals, "", reg);
  }
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();
  
  int32_t nsamples = inCrams.size();
  std::vector<std::string> v_regs;
  std::vector<int32_t> v_rids;
  std::vector<int32_t> v_poss;
  std::vector<char> v_refs;
  std::vector<char> v_alts;
  std::vector<uint8_t*> v_pls;  // stored PLs, ADs, and ODs together
  std::vector<uint8_t*> v_ads;  // stored PLs, ADs, and ODs together
  std::vector<uint8_t*> v_ods;  // stored PLs, ADs, and ODs together

  std::vector<int32_t> vSex;
  std::map<std::string,int> mSex;

  if ( !sexMap.empty() ) {
    htsFile *file = hts_open(sexMap.c_str(),"r");
    if ( file == NULL ) {
      fprintf(stderr,"ERROR: Cannot open %s\n",sexMap.c_str());
      exit(1);
    }
    kstring_t *s = &file->line;
    while( hts_getline(file,'\n',s) >= 0 ) {
      std::string ss = std::string(s->s);
      size_t idx = ss.find_first_of("\t ");
      if ( idx == std::string::npos ) {
	fprintf(stderr,"ERROR: Cannot parse line %s in %s\n",ss.c_str(), sexMap.c_str());
	exit(1);
      }
      std::string id = ss.substr(0, idx);
      int32_t sex = atoi(ss.substr(idx+1).c_str());
      
      if ( mSex.find(id) != mSex.end() ) {
	fprintf(stderr,"ERROR: Duplicate ID %s in %s\n",id.c_str(), sexMap.c_str());
	exit(1);	      
      }
      
      if ( sex == 0 ) {
	fprintf(stderr,"WARNING: Unknown sex for individual %s, assuming females\n",id.c_str());
	sex = 2;
      }
      else if ( sex > 2 ) {
	fprintf(stderr,"ERROR: Invalid sex %d for individual %s\n",sex,id.c_str());
	exit(1);
      }
      mSex[id] = sex;
    }
  }  

  notice("Started Reading site information from VCF file");
  // load all site information first
  char region[65535];
  
  while( odr.read(iv) ) {  // read marker
    bcf_unpack(iv, BCF_UN_STR);
    if ( iv->n_allele > 2 ) continue; // skip multi-allelics
    if ( !bcf_is_snp(iv) ) continue;  // focus only on SNPs

    int32_t rid = iv->rid;
    int32_t pos = iv->pos;
    char ref = iv->d.allele[0][0];
    char alt = iv->d.allele[1][0];

    sprintf(region, "%s:%d-%d", bcf_hdr_id2name(odr.hdr, rid), pos+1, pos+1);

    uint8_t* p_pls = (uint8_t*)calloc(nsamples * 3, sizeof(uint8_t));
    uint8_t* p_ads = (uint8_t*)calloc(nsamples * 2, sizeof(uint8_t));
    uint8_t* p_ods = (uint8_t*)calloc(nsamples * 1, sizeof(uint8_t));

    v_regs.push_back(region);
    v_rids.push_back(rid);
    v_poss.push_back(pos);
    v_refs.push_back(ref);
    v_alts.push_back(alt);

    v_pls.push_back(p_pls);
    v_ads.push_back(p_ads);
    v_ods.push_back(p_ods);
  }
  
  notice("Finished Reading %d site information from VCF file",(int32_t)v_poss.size());

  //#ifdef _OPENMP
  //omp_set_num_threads(threads);
  //#else
  //threads = 1;
  //#endif

  std::vector<std::string> v_sms(nsamples);
  
  // Read BAM files for each sample
  //#pragma omp parallel for schedule(dynamic, 1)
  for(int32_t i=0; i < nsamples; ++i) {  // read each CRAM files in parallel
    samFile* in = NULL;
    bam_hdr_t* hdr;
    char base, qual;
    int32_t rpos;
    
    if ( ( in = sam_open(inCrams[i].c_str(), "r") ) == 0 ) 
      error("[E:%s:%d %s] Cannot open SAM/BAM/CRAM file %s",__FILE__,__LINE__,__FUNCTION__,inCrams[i].c_str());
    
    if ( ( hdr = sam_hdr_read(in) ) == 0 )
      error("[E:%s:%d %s] Cannot open header from %s\n",__FILE__,__LINE__,__FUNCTION__,inCrams[i].c_str());

    v_sms[i] = bam_hdr_get_sample_name(hdr);
    if ( v_sms[i].empty() ) {
      size_t islash = inCrams[i].find_last_of('/');
      if ( islash == std::string::npos ) 
	v_sms[i] = inCrams[i];
      else
	v_sms[i] = inCrams[i].substr(islash+1);
      
      notice("Using %s for sample name instead", v_sms[i].c_str());
    }

    notice("Processing reg = %s, i=%d, SM=%s, nsamples %d",reg.c_str(), i,v_sms[i].c_str(),nsamples);
    
    bam1_t *b = bam_init1();

    kstring_t readseq = {0,0,0};
    kstring_t readqual = {0,0,0};	
    
    hts_idx_t *idx = sam_index_load(in, inCrams[i].c_str());
    if ( idx == NULL )
      error("[E:%s:%d %s] Cannot load index file for %s",__FILE__,__LINE__,__FUNCTION__,inCrams[i].c_str());
    int32_t numReads = 0;

    // read sam
    for(int32_t j=0; j < (int32_t)v_poss.size(); ++j) {  // process each marker separately
      //if ( j % 1000 == 0 )
      //notice("i=%d, j=%d, pos = %d", i, j, v_poss[j]);
      double p[3] = {1,1,1};
      int32_t ads[3] = {0,0,0};
      double pm, pe, sump;
      
      hts_itr_t* itr = bam_itr_querys(idx, hdr, v_regs[j].c_str());
      while( sam_itr_next(in, itr, b) >= 0 ) {
	bam_get_base_and_qual_and_read_and_qual(b, (uint32_t)v_poss[j], base, qual, rpos, &readseq, &readqual);
	//free(readseq.s);
	//free(readqual.s);
	if ( qual < 34 ) qual = 34;
	if ( qual > 73 ) qual = 73;

	pm = phredConv.phred2Mat[qual-33];
	pe = phredConv.phred2Err[qual-33];	

	if ( base == v_refs[j] ) {
	  p[0] *= ( pm * (1-gl_adj) + pe * gl_adj / 3 );
	  p[1] *= ( pm / 2 + pe / 6 );
	  p[2] *= ( pm * gl_adj + pe * (1-gl_adj) / 3 );
	  ++ads[0];
	}
	else if ( base == v_alts[j] ) {
	  p[0] *= ( pm * gl_adj + pe * (1-gl_adj) / 3 );	  	  
	  p[1] *= ( pm / 2 + pe / 6 );
	  p[2] *= ( pm * (1-gl_adj) + pe * gl_adj / 3 );	  
	  ++ads[1];
	}
	else {
	  ++ads[2];
	}
	sump = p[0]+p[1]+p[2]+1e-300;
	p[0] /= sump;
	p[1] /= sump;
	p[2] /= sump;
	
	++numReads;
	//notice("%d\t%d\t%d\t%d\t%c\t%c\t%d",i,omp_get_thread_num(),omp_get_num_threads(),v_poss[j],base,qual,rpos);
      }
      sam_itr_destroy(itr);

      sump = p[0];
      if ( p[1] > sump ) sump = p[1];
      if ( p[2] > sump ) sump = p[2];
      p[0] /= sump;
      p[1] /= sump;
      p[2] /= sump;      

      v_pls[j][i * 3 + 0] = phredConv.err2Phred(p[0]);
      v_pls[j][i * 3 + 1] = phredConv.err2Phred(p[1]);
      v_pls[j][i * 3 + 2] = phredConv.err2Phred(p[2]);
      v_ads[j][i * 2 + 0] = (ads[0] > 255 ? 255 : ads[0]);
      v_ads[j][i * 2 + 1] = (ads[1] > 255 ? 255 : ads[1]);
      v_ods[j][i] = (ads[2] > 255 ? 255 : ads[2]);
    }
    hts_idx_destroy(idx);
    sam_close(in);
    bam_hdr_destroy(hdr);
    bam_destroy1(b);    
  }

  notice("Finished reading all BAM files to calculate genotype likelihoods and allele depths");

  // write BCF output files
  BCFOrderedWriter odw(out.c_str(), 0);
  bcf_hdr_transfer_contigs((const bcf_hdr_t*) odr.hdr, odw.hdr );
  //odw->set_hdr(odr.hdr);

  odr.close();  

  bcf_hdr_append(odw.hdr, "##INFO=<ID=AVGDP,Number=1,Type=Float,Description=\"Average Depth per Sample\">\n");	
  bcf_hdr_append(odw.hdr, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate Allele Counts\">\n");
  bcf_hdr_append(odw.hdr, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total Number Allele Counts\">\n");
  bcf_hdr_append(odw.hdr, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alternate Allele Frequency from Best-guess Genotypes\">\n");
  bcf_hdr_append(odw.hdr, "##INFO=<ID=GC,Number=G,Type=Integer,Description=\"Genotype Counts\">\n");
  bcf_hdr_append(odw.hdr, "##INFO=<ID=GN,Number=1,Type=Integer,Description=\"Total Number of Genotypes\">\n");
  bcf_hdr_append(odw.hdr, "##INFO=<ID=HWEAF,Number=A,Type=Float,Description=\"Genotype likelihood based Allele Frequency assuming HWE\">\n");
  bcf_hdr_append(odw.hdr, "##INFO=<ID=HWDGF,Number=G,Type=Float,Description=\"Genotype likelihood based Genotype Frequency ignoring HWE\">\n");
  bcf_hdr_append(odw.hdr, "##INFO=<ID=IBC,Number=1,Type=Float,Description=\"Inbreeding Coefficients calculated from genotype likelihoods\">\n");	
  bcf_hdr_append(odw.hdr, "##INFO=<ID=HWE_SLRT,Number=1,Type=Float,Description=\"Signed LRT test statistics based Hardy Weinberg ln(Likelihood Ratio)\">\n");
  bcf_hdr_append(odw.hdr, "##INFO=<ID=ABE,Number=1,Type=Float,Description=\"Expected allele Balance towards Reference Allele on Heterozygous Sites\">\n");
  //bcf_hdr_append(odw.hdr, "##INFO=<ID=NS_NREF,Number=1,Type=Integer,Description=\"Number of samples with non-reference reads\">\n");
  bcf_hdr_append(odw.hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
  bcf_hdr_append(odw.hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n");
  bcf_hdr_append(odw.hdr, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele Depth\">\n");
  bcf_hdr_append(odw.hdr, "##FORMAT=<ID=OD,Number=1,Type=Integer,Description=\"Other Allele Depth\">\n");
  bcf_hdr_append(odw.hdr, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scale Genotype Likelihoods\">\n");

  for(int32_t i=0; i < nsamples; ++i) {
    bcf_hdr_add_sample(odw.hdr, v_sms[i].c_str());
  }

  odw.write_hdr();
  //odw->close();

  // write each variant in parallel
  //#pragma omp parallel for ordered schedule(static, 1) // OPENMP parallelization
  for(int32_t j=0; j < (int32_t)v_poss.size(); ++j) {  // process each marker separately
    bcf1_t* nv = bcf_init();
    nv->rid = v_rids[j];
    nv->pos = v_poss[j];
    nv->rlen = 1;
    nv->n_sample = nsamples;

    char tmp_allele_str[4] = {0,',',0,0};
    tmp_allele_str[0] = v_refs[j];
    tmp_allele_str[2] = v_alts[j];    

    //notice("Alleles are %s %s",tmp_d_alleles[0], tmp_d_alleles[1]);
    //bcf_update_alleles(odw->hdr, nv, tmp_d_alleles, 2);
    bcf_update_alleles_str(odw.hdr, nv, tmp_allele_str);
    //notice("Successfully updated alleles to %s %s",tmp_d_alleles[0], tmp_d_alleles[1]);    

    bcf_unpack(nv, BCF_UN_ALL);

    bool isX = ( ( xLabel.compare(0, xLabel.size(), v_regs[j]) == 0 ) && ( v_regs[j].at(xLabel.size()) == ':' ) ) ? true : false;
    // calculate the allele frequencies under HWE. When calculating allele frequencies, the sex information will be ignored
    float MLE_HWE_AF[2];
    float MLE_HWE_GF[3];
    int32_t ploidy = 2; // temporarily constant
    int32_t n = 0;

    // calculate the genotypes (diploid only)
    double gp, gp_sum, max_gp;
    int32_t best_gt;
    int32_t best_a1, best_a2;
    int32_t* pls_i;
    int32_t an = 0;
    int32_t acs[2];
    int32_t gcs[3];
    float afs[3];
    int32_t* gqs = (int32_t*) calloc( nsamples, sizeof(int32_t) );
    int32_t max_gq = 0;
    int32_t* gts = (int32_t*) calloc( nsamples * 2, sizeof(int32_t) );
    int32_t* pls = (int32_t*) calloc( nsamples * 3, sizeof(int32_t) );
    int32_t* ads = (int32_t*) calloc( nsamples * 2, sizeof(int32_t) );
    int32_t* ods = (int32_t*) calloc( nsamples, sizeof(int32_t) );        
    int32_t dp_sum = 0;    

    memset(acs, 0, sizeof(int32_t)*2);
    memset(gcs, 0, sizeof(int32_t)*3);

    for(int32_t i=0; i < nsamples; ++i) {
      pls[3*i  ] = (int32_t)v_pls[j][3*i];
      pls[3*i+1] = (int32_t)v_pls[j][3*i+1];
      pls[3*i+2] = (int32_t)v_pls[j][3*i+2];

      ads[2*i  ] = (int32_t)v_ads[j][2*i];
      ads[2*i+1] = (int32_t)v_ads[j][2*i+1];
      ods[i] = (int32_t)v_ods[j][i]; 

      dp_sum += (v_ads[j][2*i] + v_ads[j][2*i+1] + v_ods[j][i]);
    }
      
    Estimator * est = new Estimator();
    est->compute_gl_af_hwe(pls, nsamples, ploidy, 2, MLE_HWE_AF, MLE_HWE_GF,  n, 1e-20);

    int32_t adSumHet[2] = {0,0};
	  
    for(int32_t i=0; i < nsamples; ++i) {
      pls_i = &pls[ i * 3 ];
      
      if ( isX && (vSex[i] == 1) ) { // haploid
	max_gp = gp_sum = gp = ( est->lt->pl2prob(pls_i[0]) * MLE_HWE_AF[0] );
	best_gt = 0; best_a1 = 0; best_a2 = 0;
	for(size_t l=1; l < 2; ++l) {
	  gp = ( est->lt->pl2prob(pls_i[ l*(l+1)/2 + l]) * MLE_HWE_AF[l] );
	  gp_sum += gp;
	  if ( max_gp < gp ) {
	    max_gp = gp;
	    best_gt = l*(l+1)/2 + l; best_a1 = l; best_a2 = l;
	  }		
	}
      }
      else { // diploid
	max_gp = gp_sum = gp = ( est->lt->pl2prob(pls_i[0]) * MLE_HWE_AF[0] * MLE_HWE_AF[0] );
	best_gt = 0; best_a1 = 0; best_a2 = 0;
	for(size_t l=1; l < 2; ++l) {
	  for(size_t m=0; m <= l; ++m) {
	    gp = ( est->lt->pl2prob(pls_i[ l*(l+1)/2 + m]) * MLE_HWE_AF[l] * MLE_HWE_AF[m] * (l == m ? 1 : 2) );
	    gp_sum += gp;
	    if ( max_gp < gp ) {
	      max_gp = gp;
	      best_gt = l*(l+1)/2 + m; best_a1 = m; best_a2 = l;
	    }
	  }
	}
	if ( best_gt == 1 ) {
	  adSumHet[0] += ads[2*i];
	  adSumHet[1] += ads[2*i+1];
	}
      }
      
      double prob = 1.-max_gp/gp_sum;  // to calculate GQ
      if ( prob <= 3.162278e-26 )
	prob = 3.162278e-26;
      if ( prob > 1 )
	prob = 1;
      
      gqs[i] = (int32_t)est->lt->prob2pl(prob);
      
      if ( ( best_gt > 0 ) && ( max_gq < gqs[i] ) )
	max_gq = gqs[i];

      gts[2*i]   = ((best_a1 + 1) << 1);
      gts[2*i+1] = ((best_a2 + 1) << 1);	    
      an += 2;             // still use diploid representation of chrX for now.
      ++acs[best_a1];
      ++acs[best_a2];
      ++gcs[best_gt];
    }
    
    for(size_t i=0; i < 2; ++i) {
      afs[i] = acs[i]/(float)an;
    }
    
    bcf_update_format_int32(odw.hdr, nv, "GT", gts, nsamples * 2);
    bcf_update_format_int32(odw.hdr, nv, "GQ", gqs, nsamples );	  
    bcf_update_format_int32(odw.hdr, nv, "AD", ads, nsamples * 2);
    bcf_update_format_int32(odw.hdr, nv, "OD", ods, nsamples );	  
    bcf_update_format_int32(odw.hdr, nv, "PL", pls, nsamples * 3);
    
    float avgdp = (float)dp_sum / (float)nsamples;
    
    nv->qual = (float) max_gq;
    bcf_update_info_float(odw.hdr, nv, "AVGDP", &avgdp, 1);	  
    bcf_update_info_int32(odw.hdr, nv, "AC", &acs[1], 1);
    bcf_update_info_int32(odw.hdr, nv, "AN", &an, 1);
    bcf_update_info_float(odw.hdr, nv, "AF", &afs[1], 1);
    bcf_update_info_int32(odw.hdr, nv, "GC", gcs, 3);
    bcf_update_info_int32(odw.hdr, nv, "GN", &nsamples, 1);
    
    if (n) {
      float* MLE_HWE_AF_PTR = &MLE_HWE_AF[1];
      bcf_update_info_float(odw.hdr, nv, "HWEAF", MLE_HWE_AF_PTR, 1);
    }

    // calculate the allele frequencies under HWD	  
    float MLE_AF[2];
    float MLE_GF[3];
    n = 0;
    est->compute_gl_af(pls, nsamples, ploidy, 2, MLE_AF, MLE_GF,  n, 1e-20);
    if (n) {
      //float* MLE_AF_PTR = &MLE_AF[1];
      //bcf_update_info_float(odw->hdr, nv, "HWDAF", MLE_AF_PTR, n_alleles-1);
      bcf_update_info_float(odw.hdr, nv, "HWDGF", &MLE_GF, 3);
    }

    if ( isX && !mSex.empty() ) { // copy only female GLs to calculate IBC and HWE_SLP
      int32_t* p_XX_pls = (int32_t*) malloc(nsamples * 3 * sizeof(int32_t));
      int32_t i, k, l;
      for(i=0, k=0; i < nsamples; ++i) {
	if ( vSex[i] == 2 ) {
	  for(l=0; l < 3; ++l)  {
	    p_XX_pls[3 * k + l] = pls[3 * i + l];
	  }
	  ++k;
	}
      }
      
      float MLE_HWE_AF_XX[2];
      float MLE_HWE_GF_XX[3];
      float MLE_AF_XX[2];
      float MLE_GF_XX[3];
      
      // calculate allele frequencies using females
      est->compute_gl_af_hwe(p_XX_pls, j, ploidy, 2, MLE_HWE_AF_XX, MLE_HWE_GF_XX,  n, 1e-20);
      est->compute_gl_af(p_XX_pls, j, ploidy, 2, MLE_AF_XX, MLE_GF_XX,  n, 1e-20);
      
      for(i=0; i < 2; ++i) {
	if ( MLE_HWE_AF_XX[i] < 1e-6 ) MLE_HWE_AF_XX[i] = 1e-6;
	if ( MLE_AF_XX[i] < 1e-6 ) MLE_AF_XX[i] = 1e-6;	      
      }
      
      for(i=0; i < 3; ++i) {
	if ( MLE_HWE_GF_XX[i] < 1e-10 ) MLE_HWE_GF_XX[i] = 1e-10;
	if ( MLE_GF_XX[i] < 1e-10 ) MLE_GF_XX[i] = 1e-10;	      
      }	    
      
      float fic = 0;
      n = 0;
      est->compute_gl_fic(p_XX_pls, j, ploidy, MLE_HWE_AF_XX, 2, MLE_GF_XX, fic, n);
      if ( isnan(fic) ) fic = 0;	  
      if (n) {
	bcf_update_info_float(odw.hdr, nv, "IBC", &fic, 1);
      }

      float abe = (adSumHet[0] + 0.5)/(adSumHet[0] + adSumHet[1] + 1.0);
      bcf_update_info_float(odw.hdr, nv, "ABE", &abe, 1);      
      
      // calculate the LRT statistics related to HWE
      float lrts;
      //float logp;
      //int32_t df;
      n = 0;
      est->compute_hwe_lrt(p_XX_pls, j, ploidy, 2, MLE_HWE_GF_XX, MLE_GF_XX, n, lrts);
      if (n) {
	if ( lrts < 0 ) lrts = 0;
	if ( fic < 0 ) lrts = 0-lrts;
	bcf_update_info_float(odw.hdr, nv, "HWE_SLRT", &lrts, 1);
      }
      
      free(p_XX_pls);
    }
    else {
      float fic = 0;
      n = 0;
      est->compute_gl_fic(pls, nsamples, ploidy, MLE_HWE_AF, 2, MLE_GF, fic, n);
      if ( isnan(fic) ) fic = 0;
      if (n) {
	bcf_update_info_float(odw.hdr, nv, "IBC", &fic, 1);
      }
      
      // calculate the LRT statistics related to HWE
      float lrts;
      //float logp;
      //int32_t df;
      n = 0;
      est->compute_hwe_lrt(pls, nsamples, ploidy, 2, MLE_HWE_GF, MLE_GF, n, lrts);
      if (n) {
	if ( lrts < 0 ) lrts = 0;	
	if ( fic < 0 ) lrts = 0-lrts;
	bcf_update_info_float(odw.hdr, nv, "HWE_SLRT", &lrts, 1);
      }

      float abe = (adSumHet[0] + 0.5)/(adSumHet[0] + adSumHet[1] + 1.0);
      bcf_update_info_float(odw.hdr, nv, "ABE", &abe, 1);            
    }

    delete est;

    //#pragma omp ordered
    //notice("Writing variant j=%d",j);
    odw.write(nv);
    bcf_destroy(nv);
    free(gts);
    free(gqs);
    free(pls);
    free(ads);
    free(ods);
    delete est;
  }

  odw.close();
  
  return 0;
}

