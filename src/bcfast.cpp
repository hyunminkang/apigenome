//#include <iostream>
//#include <string>
//#include <map>
//#include <vector>
//#include <ctime>
//#include <cstdlib>
#include <getopt.h>

#include "params.h"
#include "Error.h"
#include "htslib/vcf.h"
#include "genome_interval.h"
#include "bcf_ordered_reader.h"

#define WINDOW_SIZE 65536

class pArgs {
public:
  // VCF-related string arguments
  std::string bcf;    // input VCF or VCF file
  std::string region;
  std::string rule;
  std::string field;
  std::string scoref;

  // Other input files
  std::string indf;
  std::string itvf;
  std::string bedf;
  std::string genef;
  std::string ref;

  std::string outf;
  int unit;
  bool verbose;
  bool ignoreFilter;
  bool ignoreMissing;
  bool includeMultiAllelic;

  int minAC;
  int minMAC;
  int maxAC;
  double minMAF;
  double maxMAF;
  double minCallRate;

  bool genoFlag;
  bool acFlag;
  bool anFlag;
  bool aldFlag;
  bool tstvFlag;

  bool sepchr;

  static int const DEFAULT_UNIT = 10000L;
  static double const DEFAULT_MIN = 1e-6;
  static double const DEFAULT_MAX_MAF = 1;

  pArgs() :
    unit(DEFAULT_UNIT), verbose(false), ignoreFilter(false), ignoreMissing(false), includeMultiAllelic(false), minAC(0), minMAC(0), maxAC(INT_MAX), minMAF(DEFAULT_MIN), maxMAF(DEFAULT_MAX_MAF), minCallRate(DEFAULT_MIN), genoFlag(false), acFlag(false), anFlag(false), aldFlag(false), tstvFlag(false),sepchr(false)
    {}
};

int runSummary(int argc, char** argv) {
  pArgs arg;
  
  bool biallelicOnly = false;
  bool allvar = false;
  bool snps = false;
  bool nonsnps = false;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("VCF Input Options", NULL)
    LONG_STRING_PARAM("bcf",&arg.bcf, "Input BCF or VCF file")
    LONG_STRING_PARAM("region",&arg.region, "Region to focus on in 'chr:beg-end' format")
    LONG_STRING_PARAM("indf",&arg.indf,"File containing individual IDs to focus on")
    LONG_DOUBLE_PARAM("minMAF",&arg.minMAF,"Minumum minor allele frequency cutoff")
    LONG_DOUBLE_PARAM("minCallRate",&arg.minCallRate,"Minimum call rate threshold")
    LONG_PARAM("ignore-filter",&arg.ignoreFilter,"Ignore filter and consider non-PASS variants")    
    LONG_INT_PARAM("minAC",&arg.minAC,"Minimum non-reference allele count cutoff")
    LONG_INT_PARAM("maxAC",&arg.maxAC,"Maximum non-reference allele count cutoff")
    LONG_PARAM("biallelic-only",&biallelicOnly,"Skip multi-allelic variants")    

    LONG_PARAM_GROUP("Variant Types", NULL)    
    EXCLUSIVE_PARAM("all-variants",&allvar,"Focus on all types of variants")
    EXCLUSIVE_PARAM("snps-only",&snps,"Focus on SNPs only")
    EXCLUSIVE_PARAM("exclude-snps",&nonsnps,"Focus on non-SNP variants")    

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out",&arg.outf,"Output file name")
    LONG_PARAM("verbose",&arg.verbose,"Turn on verbose mode")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( arg.bcf.empty() || arg.outf.empty()  ) {
    error("--vcf, --out are required parameters (--indf are also recommended)");
  }

  std::vector<GenomeInterval> intervals;  
  BCFOrderedReader odr(arg.bcf, intervals);
 
  bcf1_t* iv = bcf_init();
  while( odr.read(iv) ) {
    bool skip = false;
    bcf_unpack(iv, BCF_UN_ALL);

    if ( biallelicOnly && iv->n_allele > 2 ) continue;
    if ( snps || nonsnps ) {
      bool is_snp = bcf_is_snp(iv);
      if ( ( nonsnps && is_snp ) || ( snps && ( !is_snp ) ) ) continue;
    }
  }

  return 0;
}

int main(int argc, char** argv) {
  if ( argc < 2 ) {
    printf("vcfast -- Fast analytic tools for analyzing and manipulating VCF and BCF files\n");
    printf("Copyright (c) 2015 Hyun Min Kang\n");
    printf("Usage : %s [gen-hash] [options]\n",argv[0]);
    printf("\tType one of the following commands below to get detailed usage\n");
    printf("\t%s summary      [options] : summarize vcf\n",argv[0]);
  }
  else {
    std::string cmd(argv[1]);
    if ( cmd == "summary" ) {
      return runSummary(argc-1,argv+1);
    }
    else {
      error("Unrecognized command %s\n",argv[0]);
    }
  }
  return 0;
}
