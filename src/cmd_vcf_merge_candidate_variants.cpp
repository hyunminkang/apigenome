#include "cramore.h"
#include "estimator.h"
#include "htslib/kseq.h"
#include "tsv_reader.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "sex_ploidy_map.h"
#include "frequency_estimator.h"

#define MAX_SAMPLES 10

#define VT_SNP      1   //min(rlen,alen)==1 && diff==1
#define VT_INDEL    4   //diff!=0 && (rlen==1 || alen==1)

struct bcf1_comp {
  bool operator() (const bcf1_t* lhs, const bcf1_t* rhs) const {
    if ( lhs->rid == rhs->rid ) {
      if ( lhs->pos == rhs->pos ) {
	if ( lhs->rlen == rhs->rlen ) {
	  if ( lhs->n_allele == rhs->n_allele ) {
	    bcf_unpack((bcf1_t*)lhs, BCF_UN_STR);
	    bcf_unpack((bcf1_t*)rhs, BCF_UN_STR);
	    for(int32_t i=0; i < lhs->n_allele; ++i) {
	      int32_t cmp = strcmp(lhs->d.allele[i], rhs->d.allele[i]);
	      if ( cmp != 0 ) return ( cmp < 0 );
	    }
	    return false;
	  }
	  else return (lhs->n_allele < rhs->n_allele);
	}
	else return (lhs->rlen < rhs->rlen);		  
      }
      else return (lhs->pos < rhs->pos);	
    }
    else return (lhs->rid < rhs->rid);
  }
};

static std::vector<double> logfacs;
    
static double log_d(int32_t x, int32_t n, double p) {
  if ( x > n ) x = n;
  if ( ( x < 0 ) || ( n < 0 ) || ( x > n ) || ( p <= 0 ) || ( p >= 1 ) ) {
    fprintf(stderr,"[E:%s:%d %s] Cannot calculate binomial density with (%d, %d, %lf)\n", __FILE__, __LINE__, __FUNCTION__, x, n, p);
    exit(1);
  }
  for(int32_t i=(int32_t)logfacs.size(); i <= n; ++i) {
    if ( i == 0 ) logfacs.push_back(0);
    else logfacs.push_back(logfacs[i-1] + log((double)i));
  }
  return ( logfacs[n] - logfacs[x] - logfacs[n-x] + x * log(p) + (n-x) * log(1.0-p) );
}

class varMergeInfo {
public:
  int32_t nsamples;
  int32_t esum;
  float maxqual;
  
  std::vector<std::string> ids;
  int32_t ab20[20];
  int32_t dp20[20];    
  
  varMergeInfo() : nsamples(0), maxqual(0) {
    memset(ab20, 0, sizeof(int32_t)*20);
    memset(dp20, 0, sizeof(int32_t)*20);      
  }

  void add_many(int32_t ns, char* sample_ids, int32_t* _ab20, int32_t* _dp20, float qual) {
    if ( ns == 0 ) return;
    if ( qual > maxqual ) maxqual = qual;

    kstring_t s = {0,0,NULL};
    s.l = s.m = strlen(sample_ids);
    s.s = sample_ids;
    int32_t nfields = 0;
    int32_t *fields = ksplit(&s, ',', &nfields);
    if ( nsamples + ns <= MAX_SAMPLES ) {
      if ( ns < nfields )
	error("[E:%s:%d:%s] ns = %d < nfields = %d", __FILE__, __LINE__, __FUNCTION__, ns, nfields);
      for(int32_t i=0; i < nfields; ++i) {
	ids.push_back(&s.s[fields[i]]);
      }
    }
    else {
      // need to randomly sample, proportionally with
      // nsamples : ids.size() = ns : nfields
      // x = (nsamples+ns) choose MAX_SAMPLES
      std::vector<int32_t> r(nsamples+ns);
      std::fill(r.begin(), r.end(), 0);
      std::fill(r.begin() + nsamples, r.end(), 1);
      for(int32_t i=0; i < MAX_SAMPLES; ++i) {
	int32_t idx = i+(int)floor((rand()+0.5)/(RAND_MAX+1.0)*(nsamples+ns-i));
	int32_t tmp = r[idx];
	r[idx] = r[i];
	r[i] = tmp;
      }
      // shuffle ids
      std::vector<int32_t> i1(ids.size());
      for(int32_t i=0; i < (int32_t)ids.size(); ++i) i1[i] = i;
      for(int32_t i=0; i < (int32_t)ids.size(); ++i) {
	int32_t idx = i+(int32_t)floor((rand()+0.5)/(RAND_MAX+1.0)*(ids.size()-i));
	int32_t tmp = i1[idx];
	i1[idx] = i1[i];
	i1[i] = tmp;
      }

      // shuffle ids
      std::vector<int32_t> i2(nfields);
      for(int32_t i=0; i < nfields; ++i) i2[i] = i;
      for(int32_t i=0; i < nfields; ++i) {
	int32_t idx = i+(int32_t)floor((rand()+0.5)/(RAND_MAX+1.0)*(nfields-i));
	int32_t tmp = i2[idx];
	i2[idx] = i2[i];
	i2[i] = tmp;
      }
      
      std::vector<std::string> new_ids;
      int32_t j=0, k = 0;
      for(int32_t i=0; i < MAX_SAMPLES; ++i) {
	if ( r[i] == 0 ) 
	  new_ids.push_back(ids[i1[j++]]);
	else
	  new_ids.push_back(&s.s[fields[i2[k++]]]);	  
      }
      ids.swap(new_ids);
    }

    for(int32_t i=0; i < 20; ++i) {
      ab20[i] += _ab20[i];
      dp20[i] += _dp20[i];
    }
    
    nsamples += ns;
  }

  void add(std::string& id, int32_t e, int32_t n, float qual) {
    if ( n == 0 ) return;
    
    if ( e > n ) e = n;
    
    if ( qual > maxqual ) maxqual = qual;

    if ( nsamples < MAX_SAMPLES ) // store first 10 individuals
      ids.push_back(id);
    else { // calculate the probability to be included
      double u = (rand()+0.5)/(RAND_MAX+1.0);
      double p = MAX_SAMPLES/(nsamples+1.0);
      if ( u < p ) {
	int idx = (int)floor((rand()+0.5)/(RAND_MAX+1.0)*MAX_SAMPLES);
	ids[idx] = id;
      }
    }
    
    ++ab20[(int32_t)(e/(n+1e-6)*20.0)];
    ++dp20[n > 99 ? 19 : n/5];
    ++nsamples;
  }
};

// merge candidate variants that allows hierarchical merging
int32_t cmdVcfMergeCandidateVariants(int32_t argc, char** argv) {
  std::vector<std::string> inVcfs;
  std::string inVcfList;
  std::string outVcf;
  std::string region;
  int32_t minQUAL = 30;
  int32_t seed = 0;
  int32_t verbose = 100;

  paramList plst;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_MULTI_STRING_PARAM("in-vcf", &inVcfs, "VCF file name to paste")
    LONG_STRING_PARAM("in-vcf-list",&inVcfList, "File containing input VCF files in each line")
    
    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out-vcf", &outVcf, "Output VCF file name")

    LONG_PARAM_GROUP("Other options", NULL)
    LONG_STRING_PARAM("region",&region,"Genomic region to focus on")
    LONG_INT_PARAM("seed",&seed,"Random seed")
    LONG_INT_PARAM("verbose",&verbose,"Verbosity parameter as the number of variants to print out interm output")    
  END_LONG_PARAMS();
  
  plst.Add(new longParams("Available Options", longParameters));
  plst.Read(argc, argv);
  plst.Status();

  // check whether input is empty
  if ( ( inVcfs.size() > 0 ) && ( !inVcfList.empty() ) ) {
    error("[E:%s:%d %s] Cannot use --vcf and --vcf-list parameters together",__FILE__,__LINE__,__PRETTY_FUNCTION__);
  }
  else if ( ( inVcfs.empty() && inVcfList.empty() ) || outVcf.empty() ) {
    error("[E:%s:%d %s] Missing required options (--vcf or --vcf-list) and --out  parameters",__FILE__,__LINE__,__PRETTY_FUNCTION__);    
  }

  if ( seed == 0 ) { // if random seed is zero, and read from a random device to set a seed
    FILE* fd = fopen("/dev/urandom","rb");
    if (fd == NULL) {
      uint32_t t = (uint32_t)std::time(NULL);
      warning("[E:%s:%s %s] Cannot open /dev/urandom. Using %u from std::time(0)",t);      
      srand(t);
    }
    else {
      if ( fread(&seed, sizeof(uint32_t), 1, fd) == 0 )
	error("[E:%s:%d %s] Cannot read a uint32_t from /dev/urandom",__FILE__,__LINE__,__PRETTY_FUNCTION__);
      notice("Using random seed %u selected from /dev/urandom", seed);      
      srand(seed);
    }
  }

  // if input is given
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
  if ( !region.empty() ) {
    parse_intervals(intervals, "", region);
  }


  // create header
  BCFOrderedWriter* odw = new BCFOrderedWriter(outVcf, 0);

  int32_t *E = (int32_t*) malloc(1*sizeof(int32_t)); //NULL; // = NULL;
  int32_t *N = (int32_t*) malloc(1*sizeof(int32_t)); // [2] = {0,0}; // *N = NULL;
  int32_t no_E = 1, no_N = 1;
  char*  samples = NULL;
  int32_t no_samples = 0;
  int32_t *ab20 = (int32_t*) malloc(20*sizeof(int32_t));
  int32_t *dp20 = (int32_t*) malloc(20*sizeof(int32_t));
  int32_t no_ab20 = 20, no_dp20 = 20;  
  
  bcf1_t* nv = bcf_init();
  bcf1_t* wv = bcf_init();
  std::map<bcf1_t*, varMergeInfo, bcf1_comp>::iterator it;
  std::string sm_id;
  
  BCFOrderedReader* odr = NULL;
  int32_t nfiles = (int32_t)inVcfs.size();

  std::map<bcf1_t*, varMergeInfo, struct bcf1_comp> variants;  

  for (int32_t i=0; i<nfiles; ++i) {
    //////////////////////
    //i/o initialization//
    //////////////////////
    odr = new BCFOrderedReader(inVcfs[i], intervals);
    if ( i == 0 ) {
      bcf_hdr_append(odw->hdr, "##fileformat=VCFv4.2");
      bcf_hdr_transfer_contigs(odr->hdr, odw->hdr);
      bcf_hdr_append(odw->hdr, "##QUAL=Maximum variant score of the alternative allele likelihood ratio: -10 * log10 [P(Non variant)/P(Variant)] amongst all individuals.");
      bcf_hdr_append(odw->hdr, "##INFO=<ID=NSAMPLES,Number=1,Type=Integer,Description=\"Number of samples.\">");
      bcf_hdr_append(odw->hdr, "##INFO=<ID=SAMPLES,Number=.,Type=String,Description=\"Samples with evidence. (up to 10 randomly selected samples)\">");
      bcf_hdr_append(odw->hdr, "##INFO=<ID=AB20,Number=20,Type=Integer,Description=\"Histogram of allele balance using E and N (0 to 100 using 20 bins)\">");
      bcf_hdr_append(odw->hdr, "##INFO=<ID=DP20,Number=20,Type=Integer,Description=\"Histogram of depth (0-100 with 20 bins)\">");		
      odw->write_hdr();
    }

    //sm_id = bcf_hdr_get_sample_name(odr->hdr,0);
    sm_id.clear();
    
    while( odr->read(nv) ) {
      if ( ( nv->pos + 1 < intervals[0].start1 ) || ( nv->pos + 1 > intervals[0].end1 ) ) continue;
      
      bcf_unpack(nv, BCF_UN_STR);
      float variant_score = bcf_get_qual(nv);
      
      if (bcf_float_is_missing(variant_score)) {
	variant_score = 0;
      }
      
      int32_t vtype, ret;
      vtype = bcf_is_snp(nv) ? VT_SNP : VT_INDEL;
      
      if ((vtype == VT_SNP && variant_score >= minQUAL) ||
	  (vtype == VT_INDEL && variant_score >= minQUAL)) {
	if (bcf_get_info_int32(odr->hdr, nv, "NSAMPLES", &N, &no_N) >= 0) {
	  if ((ret = bcf_get_info_string(odr->hdr, nv, "SAMPLES", &samples, &no_samples)) < 0) 
	    error("[E:%s:%d %s] cannot get INFO values SAMPLES from %s\n", __FILE__, __LINE__, __FUNCTION__, inVcfs[i].c_str());

	  //notice("NSAMPLES = %d", N[0]);
	  //notice("ret = %d", ret);	  
	  //notice("SAMPLES = %s", samples);

	  if (bcf_get_info_int32(odr->hdr, nv, "AB20", &ab20, &no_ab20) < 0)
	    error("[E:%s:%d %s] cannot get INFO values AB20 from %s\n", __FILE__, __LINE__, __FUNCTION__, inVcfs[i].c_str());

	  if (bcf_get_info_int32(odr->hdr, nv, "DP20", &dp20, &no_dp20) < 0)
	    error("[E:%s:%d %s] cannot get INFO values DP20 from %s\n", __FILE__, __LINE__, __FUNCTION__, inVcfs[i].c_str());
	  
	  it = variants.find(nv);
	  
	  if ( it != variants.end() ) {
	    it->second.add_many(N[0], samples, ab20, dp20, nv->qual);
	  }
	  else {
	    variants[nv].add_many(N[0], samples, ab20, dp20, nv->qual);
	    nv = bcf_init();
	  }
	}
	else {
	  if ( sm_id.empty() )
	    sm_id = bcf_hdr_get_sample_name(odr->hdr,0);

	  if (bcf_get_format_int32(odr->hdr, nv, "E", &E, &no_E) < 0 ||
	      bcf_get_format_int32(odr->hdr, nv, "N", &N, &no_N) < 0) 
	    error("[E:%s:%d %s] cannot get format values E or N from %s\n", __FILE__, __LINE__, __FUNCTION__, inVcfs[i].c_str());

	  // assuming binomial distribution, compute likelihoods
	  double p0 = log_d(E[0], N[0], 0.01);
	  double p1 = log_d(E[0], N[0], 0.50); 
	  
	  if ( p0 <= p1 ) {
	    it = variants.find(nv);
	    if ( it != variants.end() ) 
	      it->second.add(sm_id, E[0], N[0], nv->qual);
	    else {
	      variants[nv].add(sm_id, E[0], N[0], nv->qual);
	      nv = bcf_init();
	    }
	  }
	}
      }
    }
    
    odr->close();
    delete odr;
    
    if ( (i + 1) % 100 == 0 )
      fprintf(stderr,"Finished processing %d input BCF/VCF files. Current variant count is %lu\n", i+1, variants.size());
  }
  
  fprintf(stderr,"Finished processing %d input BCF/VCF files. Current variant count is %lu\n", nfiles, variants.size());

  int32_t nvariants = 0; //(int32_t) variants.size();
  
  // print all variants;
  for(it = variants.begin(); it != variants.end(); ++it) {
    bcf_clear(wv);
    if ( it->second.nsamples > 0 ) {
      wv->rid = it->first->rid;
      wv->pos = it->first->pos;
      bcf_update_alleles(odw->hdr, wv, (const char**)it->first->d.allele, it->first->n_allele);
      
      wv->qual = it->second.maxqual;
      bcf_update_info_int32(odw->hdr, wv, "NSAMPLES", &it->second.nsamples, 1);
      std::string id = it->second.nsamples > 0 ? it->second.ids[0] : ".";
      for(int32_t i=1; i < (int32_t)it->second.ids.size() ; ++i) {
	id = id + "," + it->second.ids[i];
      }
      bcf_update_info_string(odw->hdr, wv, "SAMPLES", id.c_str());
      //bcf_update_info_int32(odw->hdr, wv, "ESUM", &it->second.esum, 0);
      bcf_update_info_int32(odw->hdr, wv, "AB20", &it->second.ab20, 20);
      bcf_update_info_int32(odw->hdr, wv, "DP20", &it->second.dp20, 20);      
      odw->write(wv);
      ++nvariants;
    }
  }

  odw->close();
  delete odw;
  
  free(E);
  free(N);
  free(ab20);
  free(dp20);
  free(samples);

  notice("Finished writing %d variants to BCF/VCF file %s", nvariants, outVcf.c_str());
  return 0;
}
