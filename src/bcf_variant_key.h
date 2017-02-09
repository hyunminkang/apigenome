#ifndef __BCF_VARIANT_KEY_H
#define __BCF_VARIANT_KEY_H

#include <string>
#include <vector>

#include "hts_utils.h"

class variantKeyS {
 public:
  std::string chrom;
  int32_t pos;
  int32_t rlen;
  std::vector<std::string> alleles;

  variantKeyS(bcf_hdr_t* hdr, bcf1_t* v) : chrom(bcf_get_chrom(hdr,v)), pos(v->pos), rlen(v->rlen) {
    alleles.resize(v->n_allele);
    bcf_unpack(v,BCF_UN_STR);    
    for(int32_t i=0; i < v->n_allele; ++i) {
      alleles[i].assign(v->d.allele[i]);
    }   
  }

  bool operator<(const variantKeyS& rhs) const {
    if ( chrom == rhs.chrom ) {
      if ( pos == rhs.pos ) {
	if ( rlen == rhs.rlen ) {
	  for(int32_t i=0; i < (int32_t)alleles.size(); ++i) {
	    if ( i > (int)rhs.alleles.size() ) return false;
	    else if ( alleles[i] != rhs.alleles[i] ) return alleles[i] < rhs.alleles[i];
	  }
	  return false;
	}
	else return ( rlen < rhs.rlen );
      }
      else return ( pos < rhs.pos );
    }
    else return ( chrom < rhs.chrom );
  }
};

#endif
