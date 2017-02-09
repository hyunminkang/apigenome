#include "Error.h"
#include "PhredHelper.h"

struct _bcf_vfilter_arg {
  std::vector<std::string> required_filters; // require at least one of the
  std::string include_expr;
  std::string exclude_expr;
};

struct _bcf_gfilter_arg {
  int32_t minDP;
  int32_t minGQ;
  //int32_t minAD;

  _bcf_gfilter_arg() {
    minDP = 0;
    minGQ = 0;
  }
};

typedef struct _bcf_vfilter_arg bcf_vfilter_arg;
typedef struct _bcf_gfilter_arg bcf_gfilter_arg;

#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

class variantKeyS {
 public:
  std::string chrom;
  int32_t pos;
  int32_t rlen;
  std::vector<std::string> alleles;

  variantKeyS(bcf_hdr_t* hdr, bcf1_t* v) : chrom(bcf_get_chrom(hdr,v)), pos(v->pos), rlen(v->rlen) {
    alleles.resize(v->n_allele);
    for(int32_t i=0; i < v->n_allele; ++i) {
      alleles[i].assign(v->d.allele[i]);
    }   
 }
};

class variantKey {
 public:
  int32_t rid;
  int32_t pos;
  int32_t rlen;
  std::vector<std::string> alleles;

  variantKey(bcf1_t* v) : rid(v->rid), pos(v->pos), rlen(v->rlen) {
    alleles.resize(v->n_allele);
    for(int32_t i=0; i < v->n_allele; ++i) {
      alleles[i].assign(v->d.allele[i]);
    }
  }
};

class variantPool {
 public:
  int32_t poolSize;
  bcf_hdr_t* hdr;
  int32_t head;
  int32_t size;
  int32_t nInfos;
  int32_t nFmts;
  std::vector<variantKey> variants;
  
  std::vector<int32_t> infoIds;
  std::vector<int32_t> infoNums;
  std::vector<std::string> infoStrs;
  std::vector< std::vector<double> > infoVals;  
  
  std::vector<int32_t> fmtIds;
  std::vector<int32_t> fmtNums;
  std::vector<std::string> fmtStrs;  
  std::vector< std::vector<double> > fmtVals;  

  variantPool(int32_t _poolSize, bcf_hdr_t* _hdr, std::vector<std::string>& infoFields, std::vector<std::string>& fmtFields) : poolSize(_poolSize), hdr(_hdr), head(0), size(0) {
    nInfos = (int32_t)infoFields.size();
    nFmts = (int32_t)fmtFields.size();

    // extract INFO fields
    for(int32_t i=0; i < nInfos; ++i) {
      uint32_t icolon = infoFields[i].find(':');
      std::string field;
      if ( icolon == std::string::npos ) {
	field = infoFields[i];
      }
      else {
	field = infoFields[i].substr(0, icolon);
      }
      int32_t id = bcf_hdr_id2int(hdr, BCF_DT_ID, field.c_str());
      if ( id < 0 )
	error("Cannot find INFO column ID %s from VCF file",field.c_str());
    }
  }
};

class InVcfArgs {
 public:
  std::string inVcf;                       // input file name
  std::string regionBed;                   // input region as BED files
  std::string regionStr;                   // input region as string
  std::vector<std::string> fmtFields;      // FORMAT fields to extract from. Example is GT, EC, DP, AD:1, AD:2, AD:SUM, AD:AVG, AD:NUM
  std::vector<std::string> infoFields;     // INFO fields to extract from. Example includes AC:1, AC:SUM, AN, ...
  bcf_vfilter_arg vfilt;
  bcf_gfilter_arg gfilt;
};
