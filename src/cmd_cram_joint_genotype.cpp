#include "cramore.h"
#include "bcf_filtered_reader.h"
#include "bcf_ordered_writer.h"
#include "frequency_estimator.h"
#include <map>
#include <string>
#include <ctime>
#include "Eigen/Dense"

//typedef std::map<std::string,double*>::iterator itU_t;

int32_t cmdCramJointGenotype(int32_t argc, char** argv) {
  std::string sampleIndex;
  std::string siteVcf;
  std::string outVcf;  
  std::string refFasta;
  double maxMemoryGB = 2.0;
  std::string tmpPrefix;
  std::string region;
  std::string sexMapFile;
  int32_t maxBQ = 40;
  int32_t minContam = 0.01;
  int32_t xStart = 2699520;
  int32_t xStop = 154931044;
  std::string xLabel("X");
  std::string yLabel("Y");
  std::string mtLabel("MT");  
  int32_t seed = 0;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Files", NULL)
    LONG_STRING_PARAM("index",&sampleIndex, "Sample index file containing [SAMPLE_ID] [CRAM_PATH] [CONTAM] [PC1] [PC2] ...")
    LONG_STRING_PARAM("site",&siteVcf, "Input VCF/BCF file containing sites to genotypes")
    LONG_STRING_PARAM("out",&outVcf, "Output VCF/BCF file")    
    LONG_STRING_PARAM("ref",&refFasta, "Reference FASTA file name (specify only when chunking is used)")    

    LONG_PARAM_GROUP("Options to specify when chunking is used", NULL)
    LONG_DOUBLE_PARAM("mem-gb",&maxMemoryGB, "Maximum (approximate) memory to use in GB")
    LONG_STRING_PARAM("tmp",&tmpPrefix, "Prefix of temporary file (default: --out parameter)")
    LONG_STRING_PARAM("region",&region, "Region to select")
    LONG_STRING_PARAM("sex-map",&sexMapFile, "File containing sex map for inidividuals (for non-autosomes)")                
    LONG_INT_PARAM("max-bq",&maxBQ, "Maximum base quality")    
    LONG_INT_PARAM("unit",&bfr.unit, "Chunking unit in bp (specify only with --ref together")
    LONG_STRING_PARAM("interval",&bfr.interval_file_name, "Interval file name used for chunking (specify only when chunking is used without --ref")
    LONG_STRING_PARAM("region",&bfr.target_region, "Target region to focus on")    

    LONG_PARAM_GROUP("Samples to focus on",NULL)
    LONG_STRING_PARAM("sm",&smID, "Sample ID to infer ancestry")    
    LONG_STRING_PARAM("sm-list",&smList,"List of sample IDs to infer ancestries")

    LONG_PARAM_GROUP("Parameters for sex chromosomes", NULL)
    LONG_STRING_PARAM("sex-map", &bfr.sexMap, "Sex map file in PED format or tsv file with [ID,SEX in X ploidy]")        
    LONG_STRING_PARAM("x-label", &bfr.xLabel, "Label for X chromosome")
    LONG_STRING_PARAM("y-label", &bfr.yLabel, "Label for Y chromosome")
    LONG_STRING_PARAM("mt-label", &bfr.mtLabel, "Label for MT chromosome")
    LONG_INT_PARAM("x-start", &bfr.xStart, "Start coordinate of non-PAR X region")
    LONG_INT_PARAM("x-stop", &bfr.xStop, "Stop coordinate of non-PAR X region")

    LONG_PARAM_GROUP("Output Files", NULL)
    LONG_STRING_PARAM("out",&outVcf, "Output VCF file to write with ISHWEZ and ISIBC statistics and IF format field")
  END_LONG_PARAMS();
  
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  notice("Analysis Started");    
  
  // sanity check of input arguments
  if ( outVcf.empty() || evecFile.empty() ) {
    error("[E:%s:%d %s] --evec, --out, are required parameters",__FILE__,__LINE__,__PRETTY_FUNCTION__);
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
      error("[E:%s:%d %s] observed %d < %d+1 columns in the file",__FILE__,__LINE__,__PRETTY_FUNCTION__, ncols, numPC);

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
  if ( ( !smID.empty() ) && ( !smList.empty() ) )
    error("[E:%s:%d %s] --sm and --sm-list cannot be used together",__FILE__,__LINE__,__PRETTY_FUNCTION__);
  if ( !smID.empty() )
    bfr.add_specified_sample(smID.c_str());
  else if ( !smList.empty() ) {
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
  Eigen::MatrixXd eV = Eigen::MatrixXd::Zero(ns, numPC);
  std::map<std::string, double*>::iterator it;  
  for(int32_t i=0; i < ns; ++i) {
    std::string sm = bfr.get_sample_id_at(i);
    it = sm2evecs.find(sm);
    if ( it == sm2evecs.end() )
      error("[E:%s:%d %s] Cannot find sample ID %s", __FILE__, __LINE__, __PRETTY_FUNCTION__, sm.c_str());
    for(int32_t j=0; j < numPC; ++j)
      eV(i,j) = it->second[j];
  }

  //std::vector< std::vector<double> > probs; // nsample * (3 * nvar) matrix
  //probs.resize(ns);

  // read genotype likelihoods
  //double* optLoadings = new double[numPC+1];

  BCFOrderedWriter odw(outVcf.c_str(), 0);
  odw.set_hdr(bfr.cdr.hdr);
  frequency_estimator freqest(&eV);
  freqest.set_hdr(odw.hdr);
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
    freqest.set_variant(nv, bfr.ploidies);

    //freqest.estimate_isaf_simplex();
    freqest.estimate_isaf_em();    
    freqest.update_variant();
    odw.write(nv);
    bcf_destroy(nv);
       
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

