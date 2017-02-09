#include "cramore.h"

#include <iostream>

#include "bcf_filtered_reader.h"

#include "Eigen/Dense"
#include "Eigen/Core"
#include "Eigen/SVD"

using namespace Eigen;

int32_t cmdVcfSVD(int32_t argc, char** argv) {
  BCFFilteredReader bfr;
  std::string outPrefix;
  int32_t maxPC = 20;
  int32_t threads = 1;

  bfr.vfilt.minMAF = 0.01;
  bfr.vfilt.minCallRate = 0.95;
  bfr.vfilt.maxAlleles = 2;
  bfr.verbose = 10000;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input VCFs", NULL)
    LONG_STRING_PARAM("in",&bfr.bcf_file_name, "Input VCF/BCF files")
    LONG_STRING_PARAM("region",&bfr.target_region, "Region to target")    

    LONG_PARAM_GROUP("Variant Filtering Options", NULL)    
    LONG_MULTI_STRING_PARAM("apply-filter",&bfr.vfilt.required_filters, "Require at least one of the listed FILTER strings")
    LONG_STRING_PARAM("include-expr",&bfr.vfilt.include_expr, "Include sites for which expression is true")
    LONG_STRING_PARAM("exclude-expr",&bfr.vfilt.exclude_expr, "Exclude sites for which expression is true")
    LONG_DOUBLE_PARAM("min-maf", &bfr.vfilt.minMAF, "Minimum minor allele frequency to filter")
    LONG_DOUBLE_PARAM("min-callrate", &bfr.vfilt.minCallRate, "Minimum call rate to filter")    

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_INT_PARAM("max-pc",&maxPC, "Maximum number of eigenvalues/eigenvectors to store")        
    LONG_STRING_PARAM("out", &outPrefix, "Output file prefix")
    LONG_INT_PARAM("verbose",&bfr.verbose, "Verbosity parameters")
    LONG_INT_PARAM("threads",&threads, "Number of threads")     
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();
  
  // sanity check of input arguments
  if ( bfr.bcf_file_name.empty() || outPrefix.empty() ) {
    error("[E:%s:%d %s] --in, --out are required parameters",__FILE__,__LINE__,__PRETTY_FUNCTION__);
  }

  bfr.init_params();

  if ( threads > 1 )
    Eigen::setNbThreads(threads);

  notice("Initializing parameters");

  std::vector<double*> v_gts;
  std::vector<double> v_afs;
  
  int32_t nsamples = bcf_hdr_nsamples(bfr.cdr.hdr);

  std::vector<std::string> varIDs;
  std::vector<std::string> smIDs;

  for(int32_t i=0; i < nsamples; ++i) {
    smIDs.push_back(bcf_hdr_sample_id(bfr.cdr.hdr, i));
  }

  notice("Starting to read variants");  

  while( bfr.read() ) {
    //notice("foo");
   
    if ( bfr.n_gts == 0 ) bfr.parse_genotypes();
    if ( bfr.n_gts != nsamples * 2 )
      error("[E:%s:%d %s] Non-diploid? %d haplotypes observed in %d samples", __FILE__,__LINE__,__PRETTY_FUNCTION__,bfr.n_gts,nsamples);
    double* gts = (double*)malloc(sizeof(double) * nsamples);
    double af = (double)(bfr.an-bfr.acs[0])/(double)bfr.an;
    for(int32_t i=0; i < nsamples; ++i) {
      if ( ( bfr.gts[2*i] == bcf_gt_missing ) || ( bfr.gts[2*i+1] == bcf_gt_missing ) ) {
	gts[i] = 0.0;
      }
      else {
	gts[i] = (bcf_gt_allele(bfr.gts[2*i]) > 0 ? 1 : 0)  + (bcf_gt_allele(bfr.gts[2*i+1]) > 0 ? 1 : 0) - af*2.0;
      }
    }
    v_afs.push_back(af);
    varIDs.push_back(bfr.get_var_ID());
    v_gts.push_back(gts);

    //if ( bfr.nRead > 10000 ) break;
  }

  int32_t nvars = (int32_t)v_gts.size();

  notice("[%s] Initializing %d by %d matrix...", __PRETTY_FUNCTION__, nvars, nsamples);

  //MatrixXd G(nvars,nsamples);
  MatrixXd G(nsamples,nvars); // use transposed from for speedy access
  for(Index i=0; i < nvars; ++i) {
    for(Index j=0; j < nsamples; ++j) {
      G(j,i) = v_gts[i][j];
    }
    free(v_gts[i]);
  }

  notice("[%s] Finishing initializtion and starting SVD..", __PRETTY_FUNCTION__);  
  //JacobiSVD<MatrixXd> svd(G, ComputeThinU | ComputeThinV);
  Eigen::BDCSVD<Eigen::MatrixXd> svd(G,Eigen::ComputeThinU | Eigen::ComputeThinV);
  notice("[%s] Finishing calculation of SVD..", __PRETTY_FUNCTION__);

  const MatrixXd& svdU = svd.matrixV();
  const VectorXd& svdD = svd.singularValues();  

  notice("[%s] Writing the UD matrix", __PRETTY_FUNCTION__);  
  htsFile* wf = hts_open((outPrefix+".fUD").c_str(), "w");
  for(int32_t i=0; i < nvars; ++i) {
    hprintf(wf,"%s\t%.5lg", varIDs[i].c_str(), v_afs[i]);
    for(int32_t j=0; j < maxPC; ++j) {
      hprintf(wf,"\t%.5lg", svdU(i,j)*svdD[j]);
    }
    hprintf(wf,"\n");
  }
  hts_close(wf);

  const MatrixXd& svdV = svd.matrixU();  
  notice("[%s] Writing the V matrix", __PRETTY_FUNCTION__);
  wf = hts_open((outPrefix+".V").c_str(), "w");
  for(int32_t i=0; i < nsamples; ++i) {
    hprintf(wf,"%s", smIDs[i].c_str());
    for(int32_t j=0; j < maxPC; ++j) {
      hprintf(wf,"\t%.5lg", svdV(i,j));
    }
    hprintf(wf,"\n");
  }
  hts_close(wf);

  notice("Analysis completed");

  return 0;
}

