#include "cramore.h"
#include "bcf_filtered_reader.h"
#include "ancestry_estimator.h"
#include <map>
#include <string>
#include <ctime>

//typedef std::map<std::string,double*>::iterator itU_t;

int32_t cmdVcfInferAncestry(int32_t argc, char** argv) {
  BCFFilteredReader bfr;  
  std::string svdPrefix;
  std::string outPrefix;
  std::string smID;
  std::string smList;
  int32_t numPC = 4;
  int32_t seed = 0;

  bfr.vfilt.maxAlleles = 2;
  bfr.verbose = 10000;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Files", NULL)
    LONG_STRING_PARAM("svd",&svdPrefix, "Prefix of SVD files (.fUD.gz, .V.gz)")
    LONG_STRING_PARAM("vcf",&bfr.bcf_file_name, "Input VCF/BCF file")
    LONG_DOUBLE_PARAM("thin",&bfr.vfilt.probThin, "Probability to thin the variants from BCF")
    LONG_INT_PARAM("seed",&seed, "Randome seed to set (default is to use clock)")

    LONG_PARAM_GROUP("Options to specify when chunking is used", NULL)    
    LONG_STRING_PARAM("ref",&bfr.ref_file_name, "Reference FASTA file name (specify only when chunking is used)")
    LONG_INT_PARAM("unit",&bfr.unit, "Chunking unit in bp (specify only with --ref together")
    LONG_STRING_PARAM("interval",&bfr.interval_file_name, "Interval file name used for chunking (specify only when chunking is used without --ref")
    LONG_STRING_PARAM("region",&bfr.target_region, "Target region to focus on")    

    LONG_PARAM_GROUP("Samples to focus on",NULL)
    LONG_STRING_PARAM("sm",&smID, "Sample ID to infer ancestry")    
    LONG_STRING_PARAM("sm-list",&smList,"List of sample IDs to infer ancestries")

    LONG_PARAM_GROUP("Output Files", NULL)
    LONG_STRING_PARAM("out",&outPrefix, "Output file prefix")
  END_LONG_PARAMS();
  
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  notice("Analysis Started");    
  
  // sanity check of input arguments
  if ( outPrefix.empty() || svdPrefix.empty() ) {
    error("[E:%s:%d %s] --svd, --out, are required parameters",__FILE__,__LINE__,__PRETTY_FUNCTION__);
  }

  srand(seed ? seed : std::time(NULL));

  // read marker loadings
  notice("Reading marker loadings");
  tsv_reader tsv_svd_u((svdPrefix+".fUD.gz").c_str());
  int32_t ncols = tsv_svd_u.read_line();
  if ( ncols < numPC + 2 )
    error("[E:%s:%d %s] observed %d < %d+2 columns",__FILE__,__LINE__,__PRETTY_FUNCTION__, ncols, numPC);    
  std::map< std::string, double* > var2u;
  while( ( ncols = tsv_svd_u.read_line() ) > 0 ) {
    if ( ncols < numPC + 2 )
      error("[E:%s:%d %s] observed %d < %d+2 columns in the file",__FILE__,__LINE__,__PRETTY_FUNCTION__, ncols, numPC);

    double* v = new double[numPC+1];
    //std::vector<double>& v = var2u[str_field_at(0)];
    //v.resize(numPC+1);
    for(int32_t i=0; i <= numPC; ++i) {
      v[i] = tsv_svd_u.double_field_at(i+1);
    }
    var2u[tsv_svd_u.str_field_at(0)] = v;
  }

  notice("Reading sample eigenvectors");  
  // read reference samples' eigenvectors
  tsv_reader tsv_svd_v((svdPrefix+".V.gz").c_str());
  ncols = tsv_svd_v.read_line();
  if ( ncols < numPC + 1 )
    error("[E:%s:%d %s] observed %d < %d+1 columns",__FILE__,__LINE__,__PRETTY_FUNCTION__, ncols, numPC);

  std::vector<std::string> refIDs;
  std::vector< std::vector<double> > matv;
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
  std::map<std::string,double*>::iterator itU;
  //itU_t itU;
  std::vector<double*> loadings;  // nvar * (numPC+1) matrix
  std::vector< std::vector<double> > probs; // nsample * (3 * nvar) matrix
  probs.resize(ns);

  // read genotype likelihoods
  while( bfr.read() ) {
    //notice("foo");
    bfr.parse_likelihoods();
    //notice("bar");    
    std::string& varID = bfr.get_var_ID();
    itU = var2u.find(varID);
    if ( itU != var2u.end() ) { // variant found
      //notice("goo");          
      loadings.push_back(itU->second);
      for(int32_t i=0; i < ns; ++i) {
	probs[i].push_back( phredConv.phred2Prob[bfr.get_likelihood_at(i*3)] );
	probs[i].push_back( phredConv.phred2Prob[bfr.get_likelihood_at(i*3+1)] );
	probs[i].push_back( phredConv.phred2Prob[bfr.get_likelihood_at(i*3+2)] );
      }
    }
  }

  notice("Estimating ancestry...");
  // Perform ancestry estimator
  ancestry_estimator ancest( &loadings, &probs, numPC, 0.5/(matv.size()));
  //double llk;
  double* optPC = new double[numPC];
  htsFile* wf = hts_open(outPrefix.c_str(),"w");
  hprintf(wf, "ID\tLLK1\tLLK0\tITER");
  //hprintf(wf, "TYPE\tID\tLLK1\tLLK0\tITER");
  for(int32_t i=0; i < numPC; ++i)
    hprintf(wf, "\tPC%d",i+1);
  hprintf(wf,"\n");

  /*
  for(int32_t i=0; i < (int32_t)matv.size(); ++i) {
    hprintf(wf,"REF\t%s\tNA\tNA\tNA", refIDs[i].c_str());
    for(int32_t j=0; j < numPC; ++j) {
      hprintf(wf,"\t%.5lf",matv[i][j]);
    }
    hprintf(wf,"\n");
    } */
  
  for(int32_t i=0; i < ns; ++i) {
    if ( i % 100 == 0 )
      notice("Performing ancestry estimation for %d-th individual %s", i, bfr.get_sample_id_at(i));
    ancest.optimizeLLK(i, optPC);
    //hprintf(wf,"TARGET\t%s\t%.5lf\t%.5lf\t%d", bfr.get_sample_id_at(i), ancest.llk1, ancest.llk0, ancest.iter);
    hprintf(wf,"%s\t%.5lf\t%.5lf\t%d", bfr.get_sample_id_at(i), ancest.llk1, ancest.llk0, ancest.iter);    
    for(int32_t j=0; j < numPC; ++j)
      hprintf(wf,"\t%.5lf", optPC[j]);
    hprintf(wf,"\n");
  }

  hts_close(wf);

  delete[] optPC;
  notice("Analysis Finished");

  return 0;
}

