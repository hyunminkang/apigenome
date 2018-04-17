#include "cramore.h"
#include "estimator.h"
#include "htslib/kseq.h"
#include "tsv_reader.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "sex_ploidy_map.h"
#include "frequency_estimator.h"

int32_t cmdVcfPasteCalls(int32_t argc, char** argv) {
  std::vector<std::string> inVcfs;
  std::string evecFile;
  std::string inVcfList;
  std::string out;
  std::string reg;
  int32_t xStart = 2699520;
  int32_t xStop = 154931044;
  std::string xLabel("X");
  std::string yLabel("Y");
  std::string mtLabel("MT");
  std::string sexMap;
  int32_t numPC = 4;
  bool skipTmpInfo = false;

  paramList plst;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_MULTI_STRING_PARAM("vcf", &inVcfs, "VCF file name to paste")
    LONG_STRING_PARAM("evec",&evecFile, "(REQUIRED) Name of eigenvector file, where each line contains [SAMPLE_ID] [PC1] [PC2] ..... The number of PCs could be larger than parameters specified by --num-PC")    
    LONG_STRING_PARAM("vcf-list",&inVcfList, "File containing input VCF files in each line")
    
    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file name")

    LONG_PARAM_GROUP("Other options", NULL)
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")
    LONG_STRING_PARAM("sex-map",&sexMap, "Sex map file, containing ID and sex (1 for male and 2 for female) for each individual")
    LONG_INT_PARAM("num-pc",&numPC, "Number of principal componentds to be used from the file specified by --evec ")
    LONG_PARAM("skip-tmp-info",&skipTmpInfo, "Skip writing temporary INFO field FLT20")    

    LONG_PARAM_GROUP("Sex Chromosomes",NULL)
    LONG_STRING_PARAM("xLabel", &xLabel, "Contig name for X chromosome")
    LONG_STRING_PARAM("yLabel", &yLabel, "Contig name for Y chromosome")
    LONG_STRING_PARAM("mtLabel", &mtLabel, "Contig name for MT chromosome")
    LONG_INT_PARAM("xStart", &xStart, "Start base position of non-PAR region in X chromosome")
    LONG_INT_PARAM("xStop",  &xStop,  "End base position of non-PAR region in X chromosome")    
  END_LONG_PARAMS();
  
  plst.Add(new longParams("Available Options", longParameters));
  plst.Read(argc, argv);
  plst.Status();

  if ( ( inVcfs.size() > 0 ) && ( !inVcfList.empty() ) ) {
    error("[E:%s:%d %s] Cannot use --vcf and --vcf-list parameters together");
  }
  else if ( ( inVcfs.empty() && inVcfList.empty() ) || out.empty() || ( evecFile.empty() && numPC > 0 ) ) {
    error("[E:%s:%d %s] Missing required options (--vcf or --vcf-list) and --out, --evec parameters");    
  }

  std::map<std::string, double*> sm2evecs;
  if ( numPC > 0 ) {
    notice("Reading sample eigenvectors");  
    // read reference samples' eigenvectors
    tsv_reader tsv_svd_v(evecFile.c_str());
    int32_t ncols; 
    
    while( ( ncols = tsv_svd_v.read_line() ) > 0 ) {
      if ( ncols < numPC + 1 )
	error("[E:%s:%d %s] observed %d < %d+1 columns in the file",__FILE__,__LINE__,__PRETTY_FUNCTION__, ncols, numPC);
      
      std::string smID = tsv_svd_v.str_field_at(0);
      double* v = new double[numPC];
      for(int32_t i=0; i < numPC; ++i) {
	v[i] = tsv_svd_v.double_field_at(i+1);
      }
      if ( sm2evecs.find(smID) != sm2evecs.end() )
	error("Duplicate ID %s in evecFile %s", smID.c_str(), evecFile.c_str());
      sm2evecs[smID] = v;
    }
    notice("Finished reading eigenvectors for %d individuals",sm2evecs.size());
  }
  
  if ( !inVcfList.empty() ) {
    htsFile* fp = hts_open(inVcfList.c_str(),"r");
    if ( fp == NULL )
      error("[E:%s:%d %s] Cannot open file %s for reading",__FILE__,__LINE__,__FUNCTION__, inVcfList.c_str());
  
    int32_t lstr = 0;
    int32_t* fields = NULL;
    int32_t n = 0;
    kstring_t str = {0,0,0};      
    while( ( lstr = hts_getline(fp, KS_SEP_LINE, &str) ) >= 0 ) {
      fields = ksplit(&str, 0, &n);

      if ( n != 1 )
	error("Only expecting one token in each line of %s",inVcfList.c_str());

      inVcfs.push_back(str.s + fields[0]);
    }
    hts_close(fp);
  }

  std::vector<GenomeInterval> intervals;
  if ( !reg.empty() ) {
    parse_intervals(intervals, "", reg);
  }

  BCFOrderedReader** odrs = (BCFOrderedReader**)malloc(sizeof(BCFOrderedReader*) * inVcfs.size());
  notice("Loading input VCFs...");
  for(int32_t i=0; i < (int32_t)inVcfs.size(); ++i) {
    odrs[i] = new BCFOrderedReader(inVcfs[i], intervals);
  }
  notice("Finished loading %u input VCFs files", inVcfs.size());

  // create header
  BCFOrderedWriter* odw = new BCFOrderedWriter(out);
  // copy the first header
  odw->hdr = bcf_hdr_dup(odrs[0]->hdr);
  //notice("i=%d, nsamples = %d",0, bcf_hdr_nsamples(odw->hdr));
  // and add samples for the rest
  for(int32_t i=1; i < (int32_t)inVcfs.size(); ++i) {
    int32_t nsamples = bcf_hdr_nsamples(odrs[i]->hdr);
    //notice("i=%d, nsamples = %d",i, nsamples);    
    for(int32_t j=0; j < nsamples; ++j) {
      bcf_hdr_add_sample(odw->hdr, bcf_hdr_get_sample_name(odrs[i]->hdr, j));
    }
  }
  bcf_hdr_add_sample(odw->hdr, NULL);
  bcf_hdr_sync(odw->hdr);
  
  //notice("i=x, nsamples = %d", bcf_hdr_nsamples(odw->hdr));  

  sex_ploidy_map spmap(xLabel, yLabel, mtLabel, xStart, xStop);
  spmap.load_sex_map_file(sexMap.empty() ? NULL : sexMap.c_str(), odw->hdr);    

  notice("Finished writing header of the output file %s", out.c_str());

  // check whether everyone has eigenvectors, and store them in order
  int32_t nsamples = bcf_hdr_nsamples(odw->hdr);
  //notice("nsamples = %d",nsamples);
  Eigen::MatrixXd eV = Eigen::MatrixXd::Constant(nsamples, numPC+1, 1.0);
  if ( numPC > 0 ) {
    std::map<std::string, double*>::iterator it;    
    for(int32_t i=0; i < nsamples; ++i) {
      std::string sm(bcf_hdr_get_sample_name(odw->hdr, i));
      it = sm2evecs.find(sm);
      if ( it == sm2evecs.end() ) {
	error("[E:%s:%d %s] Cannot find sample ID %s", __FILE__, __LINE__, __PRETTY_FUNCTION__, sm.c_str());  
      }
      for(int32_t j=0; j < numPC; ++j)
	eV(i,j+1) = it->second[j];    
    }
  }

  // read each variant
  bcf1_t** rvs = (bcf1_t**)malloc(sizeof(bcf1_t*) * inVcfs.size());
  for(int32_t i=0; i < (int32_t)inVcfs.size(); ++i) {
    rvs[i] = bcf_init1();
  }

  frequency_estimator freqest(&eV);
  freqest.field = "PL";
  freqest.skipIf = true;
  freqest.skipInfo = false;
  freqest.siteOnly = false;
  freqest.set_hdr(odw->hdr, odw->hdr);
  
  odw->write_hdr();

  //int32_t* gt = (int32_t*) calloc ( nsamples * 2, sizeof(int32_t) );
  int32_t* pl = (int32_t*) calloc ( nsamples * 3, sizeof(int32_t) );
  int32_t* ad = (int32_t*) calloc ( nsamples * 2, sizeof(int32_t) );
  int32_t* td = (int32_t*) calloc ( nsamples * 1, sizeof(int32_t) );
  //int32_t* gq = (int32_t*) calloc ( nsamples * 1, sizeof(int32_t) );  

  int32_t nvariants = 0;

  float* flt20 = NULL;
  float* avgdp = NULL;
  int32_t n_flt20 = 0;
  int32_t n_avgdp = 0;
  
  while( odrs[0]->read(rvs[0]) ) {
    // read a record from each VCF
    for(int32_t i=1; i < (int32_t)inVcfs.size(); ++i) {
      if ( !odrs[i]->read(rvs[i]) )
	error("Cannot read variant at i=%d", i);
      if ( ( rvs[0]->rid != rvs[i]->rid ) || ( rvs[0]->pos != rvs[i]->pos ) ) {
	error("The variant are not identical at %d:%d vs %d:%d", rvs[0]->rid, rvs[i]->rid, rvs[0]->pos, rvs[0]->pos);
      }
    }

    // check if the variant POS is located within the intervals
    if ( ( rvs[0]->pos+1 < intervals[0].start1 ) || ( rvs[0]->pos+1 > intervals[0].end1 ) )
      continue;

    // get the variant site info
    bcf1_t* nv = bcf_init();
    bcf_set_n_sample(nv, nsamples);
    bcf_set_rid(nv, rvs[0]->rid);
    bcf_set_pos1(nv, rvs[0]->pos+1);

    for(int32_t i=0; i < (int32_t)inVcfs.size(); ++i) 
      bcf_unpack(rvs[i],BCF_UN_ALL);

    kstring_t alleles = {0,0,0};
    char** tmp_alleles = bcf_get_allele(rvs[0]);
    for (size_t i=0; i< bcf_get_n_allele(rvs[0]); ++i) {
      if (i) kputc(',', &alleles);
      kputs(tmp_alleles[i], &alleles);
    }    
    bcf_update_alleles_str(odw->hdr, nv, alleles.s);

    // variables for computing INFO and FORMAT fields
    float flt20_sum[20];
    double dp_sum = 0;
    int32_t ioffset = 0;
    int32_t n_pls = 3*nsamples;
    int32_t n_ads = 2*nsamples;
    int32_t n_dps = nsamples;
    
    nv->qual = 0;

    // process each VCF separately
    for(int32_t i=0; i < (int32_t)inVcfs.size(); ++i) {
      //float flt20[21];
      //float* flt20 = NULL;
      //float avgdp;
      //float abe;
      //int32_t nflt = 0;
      //int32_t navgdp = 1;

      //notice("Reading FLT20\n");
      
      if ( bcf_get_info_float(odrs[i]->hdr, rvs[i], "FLT20", &flt20, &n_flt20) < 0 )
	error("Cannot read FLT20 field from %s at %d:%d", inVcfs[i].c_str(), rvs[i]->rid, rvs[i]->pos);

      //notice("Finished reading FLT20");

      if ( bcf_get_info_float(odrs[i]->hdr, rvs[i], "AVGDP", &avgdp, &n_avgdp) < 0 )
	error("Cannot read AVGDP field from %s at %d:%d", inVcfs[i].c_str(), rvs[i]->rid, rvs[i]->pos);
      
      dp_sum += ( (double)avgdp[0] * bcf_hdr_nsamples(odrs[i]->hdr) );

      //notice("avgdp = %f",avgdp[0]);            

      // update FLT20 INFO field
      for(int32_t j=0; j < 20; ++j) {
	if ( i == 0 ) flt20_sum[j] = flt20[j];
	else flt20_sum[j] += flt20[j];
      }

      // update QUAL
      //if ( nv->qual < rvs[i]->qual )      
      //nv->qual = rvs[i]->qual;

      // read PL, AD, DP
      void* p = (void*)&pl[ioffset*3];
      if ( bcf_get_format_int32(odrs[i]->hdr, rvs[i], "PL", &p, &n_pls) < 0 )
	error("Cannot read PL field from %s at %d:%d", inVcfs[i].c_str(), rvs[i]->rid, rvs[i]->pos);
      p = (void*)&ad[ioffset*2];
      if ( bcf_get_format_int32(odrs[i]->hdr, rvs[i], "AD", &p, &n_ads) < 0 )
	error("Cannot read AD field from %s at %d:%d", inVcfs[i].c_str(), rvs[i]->rid, rvs[i]->pos);
      p = (void*)&td[ioffset];
      if ( bcf_get_format_int32(odrs[i]->hdr, rvs[i], "DP", &p, &n_dps) < 0 )
	error("Cannot readDP field from %s at %d:%d", inVcfs[i].c_str(), rvs[i]->rid, rvs[i]->pos);
      
      ioffset += bcf_hdr_nsamples(odrs[i]->hdr);

      //free(flt20);
    }

    freqest.set_variant(nv, spmap.get_ploidies(nv), pl);
    freqest.estimate_isaf_em();
    freqest.score_test_hwe(true);
    
    // we need to add GT and GQ fields
    // update INFO fields
    avgdp[0] = (float)(dp_sum / (float)nsamples);
    bcf_update_info_float(odw->hdr, nv, "AVGDP", avgdp, 1);
    freqest.update_gt_gq(true);
    bcf_update_format_int32(odw->hdr, nv, "AD", ad, nsamples*2);
    bcf_update_format_int32(odw->hdr, nv, "DP", td, nsamples);
    bcf_update_format_int32(odw->hdr, nv, "PL", pl, nsamples*3);        
    freqest.update_variant();

    float abe = flt20_sum[16]/(flt20_sum[17]+1e-6);
    float abz = flt20_sum[18]/sqrt(flt20_sum[19]+1e-6);
    float bqz = flt20_sum[0]/sqrt(flt20_sum[1]+1e-6);
    float mqz = flt20_sum[2]/sqrt(flt20_sum[3]+1e-6);
    float cyz = flt20_sum[4]/sqrt(flt20_sum[5]+1e-6);
    float stz = flt20_sum[6]/sqrt(flt20_sum[7]+1e-6);
    float nmz = flt20_sum[8]/sqrt(flt20_sum[9]+1e-6);
    float ior = log(flt20_sum[10]/flt20_sum[11]+1e-6)/log(10.);
    float nm0 = flt20_sum[12]/(flt20_sum[13]+1e-6);
    float nm1 = flt20_sum[14]/(flt20_sum[15]+1e-6);

    bcf_update_info_float(odw->hdr, nv, "ABE", &abe, 1);
    bcf_update_info_float(odw->hdr, nv, "ABZ", &abz, 1);
    bcf_update_info_float(odw->hdr, nv, "BQZ", &bqz, 1);
    bcf_update_info_float(odw->hdr, nv, "MQZ", &mqz, 1);
    bcf_update_info_float(odw->hdr, nv, "CYZ", &cyz, 1);
    bcf_update_info_float(odw->hdr, nv, "STZ", &stz, 1);
    bcf_update_info_float(odw->hdr, nv, "NMZ", &nmz, 1);
    bcf_update_info_float(odw->hdr, nv, "IOR", &ior, 1);
    bcf_update_info_float(odw->hdr, nv, "NM0", &nm1, 1);
    bcf_update_info_float(odw->hdr, nv, "NM1", &nm0, 1);

    if ( !skipTmpInfo )
      bcf_update_info_float(odw->hdr, nv, "FLT20", flt20_sum, 20);    

    odw->write(nv);
    bcf_destroy(nv);
    ++nvariants;

    if ( nvariants % ( (int)( 10000000 / nsamples ) )  == 0 )
      notice("Writing %d variants to BCF/VCF file %s", nvariants, out.c_str());    
  }

  free(flt20);
  free(avgdp);

  odw->close();
  delete odw;

  notice("Finished writing %d variants to BCF/VCF file %s", nvariants, out.c_str());
  return 0;
}
