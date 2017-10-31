#include "cramore.h"
#include "bcf_filtered_reader.h"
#include "bcf_ordered_writer.h"
#include "frequency_estimator.h"
#include "Eigen/Dense"
#include <map>
#include <string>
#include <ctime>

//typedef std::map<std::string,double*>::iterator itU_t;

int32_t cmdVcfInferISAF(int32_t argc, char** argv) {
  BCFFilteredReader bfr;  
  std::string evecFile;
  std::string outVcf;
  //std::string smID;
  std::string smList;
  int32_t numPC = 4;
  int32_t seed = 0;
  bool skipIf   = false;
  bool skipInfo = false;
  bool siteOnly = false;
  std::string field;
  double gtError = 0.005;

  bfr.vfilt.maxAlleles = 2;
  bfr.verbose = 100;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Options", NULL)
    LONG_STRING_PARAM("evec",&evecFile, "(REQUIRED) Name of eigenvector file, where each line contains [SAMPLE_ID] [PC1] [PC2] ..... The number of PCs could be larger than parameters specified by --num-PC")
    LONG_STRING_PARAM("vcf", &bfr.bcf_file_name, "(REQUIRED) Input VCF/BCF file")
    LONG_DOUBLE_PARAM("thin", &bfr.vfilt.probThin, "Probability to randomly sample variants from BCF")
    LONG_INT_PARAM("seed",&seed, "Random seed to set (default is to use the clock time)")
    LONG_INT_PARAM("num-pc",&numPC, "Number of principal componentds to be used from the file specified by --evec ")
    LONG_STRING_PARAM("field",&field, "FORMAT field in VCF to extract the genotype likelihood or genotypes from. Only PL, GL, GT are allowed currently")
    LONG_DOUBLE_PARAM("gt-error",&gtError, "Error rates for GT field when --field GT option is used. Ignored for other fields")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out",&outVcf, "(REQUIRED) Output VCF file to write with ISHWEZ and ISIBC statistics and IF format field")
    LONG_PARAM("skip-if", &skipIf,   "Skip writing individual-specific allele frequency for each sample in output VCF/BCF")
    LONG_PARAM("skip-info", &skipInfo,   "Skip updating INFO field for each sample in output VCF/BCF")
    LONG_PARAM("site-only", &siteOnly,   "Do not write genotype information, and writes only site information (up to INFO field) in output VCF/BCF")

    LONG_PARAM_GROUP("Samples to focus on",NULL)
    //LONG_STRING_PARAM("sm",&smID, "Sample ID to subset from VCF/BCF when estimating ISAF. HWE statistics would not be meaningful in this case")
    LONG_STRING_PARAM("sm-list",&smList,"A file containg the list of sample IDs to subset")

    LONG_PARAM_GROUP("Parameters for sex chromosomes", NULL)
    LONG_STRING_PARAM("sex-map", &bfr.sexMap, "Sex map file in PED format or tsv file with [ID,SEX in X ploidy]")        
    LONG_STRING_PARAM("x-label", &bfr.xLabel, "Label for X chromosome")
    LONG_STRING_PARAM("y-label", &bfr.yLabel, "Label for Y chromosome")
    LONG_STRING_PARAM("mt-label", &bfr.mtLabel, "Label for MT chromosome")
    LONG_INT_PARAM("x-start", &bfr.xStart, "Start coordinate of non-PAR X region")
    LONG_INT_PARAM("x-stop", &bfr.xStop, "Stop coordinate of non-PAR X region")

    LONG_PARAM_GROUP("Options to specify when chunking is used", NULL)    
    LONG_STRING_PARAM("ref",&bfr.ref_file_name, "Reference FASTA file name (required only when chunking is used)")
    LONG_INT_PARAM("unit",&bfr.unit, "Chunking unit in bp (specify only with --ref together")
    LONG_STRING_PARAM("interval",&bfr.interval_file_name, "Interval file name used for chunking (specify only when chunking is used without --ref")
    LONG_STRING_PARAM("region",&bfr.target_region, "Target region to focus on")        
  END_LONG_PARAMS();
  
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  notice("Analysis Started");    
  
  // sanity check of input arguments
  if ( outVcf.empty() || evecFile.empty() || bfr.bcf_file_name.empty() ) {
    error("[E:%s:%d %s] --evec, --out, --vcf are required parameters",__FILE__,__LINE__,__PRETTY_FUNCTION__);
  }

  srand(seed ? seed : std::time(NULL));

  notice("Reading sample eigenvectors");  
  // read reference samples' eigenvectors
  tsv_reader tsv_svd_v(evecFile.c_str());
  int32_t ncols; // = tsv_svd_v.read_line();
  //if ( ncols < numPC )
  //error("[E:%s:%d %s] observed %d < %d+1 columns",__FILE__,__LINE__,__PRETTY_FUNCTION__, ncols, numPC);

  std::map<std::string, double*> sm2evecs;
  while( ( ncols = tsv_svd_v.read_line() ) > 0 ) {
    if ( ncols < numPC + 1 )
      error("[E:%s:%d %s] observed %d < %d+1 columns in the file %s line number %d",__FILE__,__LINE__,__PRETTY_FUNCTION__, ncols, numPC, evecFile.c_str(), tsv_svd_v.nlines);

    std::string smID = tsv_svd_v.str_field_at(0);
    double* v = new double[numPC];
    for(int32_t i=0; i < numPC; ++i) {
      v[i] = tsv_svd_v.double_field_at(i+1);
    }
    sm2evecs[smID] = v;
  }

  notice("Identifying sample columns to extract..");    
  // identify sample columns to extract
  std::vector<int32_t> isamples;
  if ( !smList.empty() ) {
    tsv_reader tsv_sm(smList.c_str());
    while ( ( ncols = tsv_sm.read_line() ) > 0 ) {
      bfr.add_specified_sample(tsv_sm.str_field_at(0));
    }
  }

  notice("Reading in BCFs...");      
  // initialize BCF reader
  bfr.init_params();

  //notice("Finished initizliaing BCF");        

  int32_t ns = bfr.get_nsamples();
  Eigen::MatrixXd eV = Eigen::MatrixXd::Constant(ns, numPC+1, 1.0);
  std::map<std::string, double*>::iterator it;  
  for(int32_t i=0; i < ns; ++i) {
    std::string sm = bfr.get_sample_id_at(i);
    it = sm2evecs.find(sm);
    if ( it == sm2evecs.end() )
      error("[E:%s:%d %s] Cannot find sample ID %s", __FILE__, __LINE__, __PRETTY_FUNCTION__, sm.c_str());
    for(int32_t j=0; j < numPC; ++j)
      eV(i,j+1) = it->second[j];
  }

  //std::vector< std::vector<double> > probs; // nsample * (3 * nvar) matrix
  //probs.resize(ns);

  // read genotype likelihoods
  //double* optLoadings = new double[numPC+1];

  BCFOrderedWriter odw(outVcf.c_str(), 0);
  if ( siteOnly ) {
    bcf_hdr_t* hnull = bcf_hdr_subset(bfr.cdr.hdr, 0, 0, 0);
    bcf_hdr_remove(hnull, BCF_HL_FMT, NULL);
    odw.set_hdr(hnull);
  }
  else {
    odw.set_hdr(bfr.cdr.hdr);
  }

  frequency_estimator freqest(&eV);
  // assign command arguments
  freqest.field = field;
  freqest.gtError = gtError;
  freqest.skipIf = skipIf;
  freqest.skipInfo = skipInfo;
  freqest.siteOnly = siteOnly;
  
  freqest.set_hdr(bfr.cdr.hdr, odw.hdr);
  odw.write_hdr();

  /*
  htsFile* wf = hts_open(outPrefix.c_str(), "w");
  hprintf(wf, "VARIANT\tLLK1\tLLK0\tLLKDIFF\tITER\tOLDZ\tNEWZ\tOLDAF\tNEWAF\tMINAF\tMAXAF\n");
  */

  
  while( bfr.read() ) {
    //notice("foo");
    //bfr.parse_likelihoods();

    //notice("bar");

    bcf1_t* nv = bcf_dup(bfr.cursor());
    freqest.set_variant(nv, bfr.ploidies); //, NULL, bfr.sm_icols.empty() ? NULL : &bfr.sm_icols);

    //notice("bar");    

    //freqest.estimate_isaf_simplex();
    freqest.estimate_isaf_em();
    freqest.score_test_hwe(true);    
    freqest.update_variant();

    //notice("car %d %d",bcf_hdr_nsamples(odw.hdr),nv->n_sample);    

    odw.write(nv);

    //notice("dar");
    
    bcf_destroy(nv);

    //notice("far");        

    //error("%d",bcf_hdr_nsamples(odw.hdr));    

       
    //frequency_estimator frqest( &evecs, &bfr, numPC, 0.5/ns );

    //double oldaf = (bfr.an - bfr.acs[0] + 0.5)/(bfr.an+1.0);
 

    //std::string msg;
    //catprintf(msg, "%s : (%.5lg) -> (", bfr.get_var_ID().c_str(), (bfr.an - bfr.acs[0] + 0.5)/(bfr.an+1.0)*2);
    
    //double llk = frqest.optimizeLLK(optLoadings);
    //hprintf(wf,"%s\t%.3lf\t%.3lf\t%.3lf\t%d\t%.3lf\t%.3lf\t%.3lg\t%.3lg\t%.3lg\t%.3lg\n",bfr.get_var_ID().c_str(), frqest.llk1, frqest.llk0, frqest.llk1-frqest.llk0, frqest.iter, frqest.hwe0z, frqest.hwe1z, oldaf, frqest.meanISAF, frqest.minISAF, frqest.maxISAF);

    //notice("goo");  
    
    //for(int32_t i=0; i < numPC+1; ++i) {
    //  catprintf(msg, "%s%.5lg", i == 0 ? "" : ", ", optLoadings[i]);
    //}
    //notice("%s, LLK1=%.5lf, LLK0=%.5lf, niter=%d, HWE1=%.5lf, HWE0=%.5lf, maxISAF=%.5lg, minISAF=%.5lg)", msg.c_str(), frqest.llk1, frqest.llk0, frqest.iter, frqest.hwe0z, frqest.hwe1z, frqest.maxISAF, frqest.minISAF);
  }
  //hts_close(wf);
  odw.close();


  //delete[] optLoadings;
  notice("Analysis Finished");

  return 0;
}
