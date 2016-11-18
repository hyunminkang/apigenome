#ifndef __APIGENOME_FASTA_H
#define __APIGENOME_FASTA_H

#include "htslib/hts.h"

class fastaMask {
public:
  std::vector<std::string> seqnames;
  std::vector<int32_t> seqlens;
  std::vector<char*> seqs;
  std::map<std::string,int32_t> seq2tid;
  faidx_t* fai;

  fastaMask(const char* fastaFile);
  ~fastaMask() {
    for(int32_t i=0; i < (int32_t)seqnames.size(); ++i) {
      free(seqs[i]);
    }
    fai_destroy(fai);
    
  }
  void maskRegion(const char* chrom, int32_t beg1, int32_t end0, char c = '\0');
  void maskRegion(int32_t tid, int32_t beg1, int32_t end0, char c = '\0');
  void maskVcf(const char* vcfFile);

  char getMask(const char* chrom, int32_t pos1);
  char getMask(int32_t tid, int32_t pos1);
};

#endif 
